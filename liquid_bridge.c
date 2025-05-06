// Author: Radu Cimpeanu
// Date: 3/11/2024

#define FILTERED
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "vof.h"                     // solve using VoF method
#include "contact-embed.h"           // IMPORTANT addition to enable contact angles on embedded solids
#include "fractions.h"               // initially define fraction of fluid phases
#include "view.h"                    // need to make the animations
#include "tag.h"                     // helps track droplet properties
#include "draw.h"

#define MAX_LEVEL   8                // max grid resolution (2^MAX_LEVEL)

// dimensional quantities (fluids at 20 degrees C)
#define rhoLiquid   961.0            // liquid phase density (kg/m^3)
#define rhoGas      1.21             // gas phase density (kg/m^3)

#define muLiquid    1.057            // liquid dynamic viscosity (kg/ms)
#define muGas       1.81e-2          // gas dynamic viscosity (kg/ms)

#define sig         0.039            // surface tension coefficient (N/m)

#define g_accel     9.81             // gravitational acceleration (m/s^2)

#define l_ref       1.0e-3           // cylinder radius (m)


// IMPORTANT Several non-dimensionalisation options included below, currently gravity-based
// #define u_ref       sqrt(sig/(rhoLiquid*l_ref))      // reference velocity (surface tension)
// #define u_ref       muLiquid/(rhoLiquid*l_ref)       // reference velocity (viscosity)
#define u_ref          sqrt(g_accel*l_ref)              // reference velocity (gravity)

// relevant dimensionless quantities

#define rho_ratio   rhoGas/rhoLiquid                         // density ratio
#define mu_ratio    muGas/muLiquid                           // viscosity ratio 

#define invRe       muLiquid/(rhoLiquid*u_ref*l_ref)         // inverse Reynolds number 
#define invWe       sig/(rhoLiquid*pow(u_ref,2.0)*l_ref)     // inverse squared Weber number 
#define Fr          (g_accel*l_ref/pow(u_ref,2.0))           // inverse squared Froude number

#define domainSize 3.6 /* Physical size of the box */

#define tEnd 100.001 // Timescale

vector h[];
face vector av[];

FILE * fp_stats;
FILE * fp_vol;
FILE * fp_interface;

// Top: outflow
u.n[top]  = neumann(0.);
p[top]   = dirichlet(0.);
pf[top]  = dirichlet(0.);

// Right: outflow
u.n[right]  = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

// Bottom: outflow
u.n[bottom]  = neumann(0.);
p[bottom]   = dirichlet(0.);
pf[bottom]  = dirichlet(0.);

// Left: symmetry
u.t[left] = neumann(0.);
u.n[left] = dirichlet(0.);

int maxLevel; // = 11;
double theta0; // =30.0;

int main(int argc, char * argv[])
{

  theta0 = atof(argv[1]);    // prescribed simulation end time
  maxLevel = atof(argv[2]);  // prescribed maximum resolution level
  
  init_grid (1 << maxLevel);
  size(domainSize);                     
  origin (0.0, 0.0);
  
  //print dimensionless numbers
  fprintf(stdout, "Dimensionless viscosity = %0.9f \n", invRe); fflush(stdout);
  fprintf(stdout, "Dimensionless surface tension coefficient = %0.9f \n", invWe); fflush(stdout);
  fprintf(stdout, "Dimensionless gravity = %0.9f \n", Fr); fflush(stdout);
  fprintf(stdout, "Density ratio = %0.9f \n", rho_ratio); fflush(stdout);
  fprintf(stdout, "Viscosity ratio = %0.9f \n", mu_ratio); fflush(stdout);
  fprintf(stdout, "Final time (dimensional) = %0.9f \n", tEnd*l_ref/u_ref); fflush(stdout);
  
  rho1 = 1.;
  rho2 = 1./rho_ratio;

  mu1 = invRe;
  mu2 = mu_ratio*mu1;

  a = av;

  f.sigma = invWe;
  
  // IMPORTANT Contact angle imposition
  const scalar c[] = theta0*pi/180.;
  contact_angle = c;
  
  // Pointer of the file for interface extrema
  {
    char name[200];
    sprintf(name, "loginterface.dat");
    fp_interface = fopen(name, "w");
  }
  
  // Pointer of the file to save stats
  {
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");
  }

  // Pointer of the file to save mass conservation info
  {
    char name[200];
    sprintf(name, "logvol.dat");
    fp_vol = fopen(name, "w");
  }

  DT = 1.0e-3;
  NITERMIN = 1; // default 1
  NITERMAX = 200; // default 100
  TOLERANCE = 1e-4; // default 1e-3
  
  run();

  fclose(fp_interface);
  fclose(fp_stats);
  fclose(fp_vol);
}

event acceleration (i++) {
  foreach_face(x)  
    av.x[] += 0.0;
  foreach_face(y)  
    av.y[] -= Fr;
}


event init (t = 0.)
{
  // IMPORTANT Embedded solid creation
  vertex scalar phi[];
  foreach_vertex() 
    phi[] = (sq(x - 1.4968) + sq(y - (domainSize/2.0)) - sq(1.0));
  boundary ({phi});
  fractions (phi, cs, fs);
  
  // IMPORTANT Initial interface shape given by ellipse
  fraction (f, - (sq(x - 0.0)/sq(2.2) + sq(y - (domainSize/2.0))/sq(0.784) - sq(1.0)));
  
  foreach() {
    u.x[] = 0.;
    u.y[] = 0.;
  }
}

event gfsview (t = 0.0; t += 25.0; t <= tEnd) {
    char name_gfs[200];
    sprintf(name_gfs,"Slices/Bridge-%0.1f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);
}
  
scalar omega[], omega_out[];

event movies (t += 0.1)
{
  char timestring[100];
  
  vorticity (u, omega);
  
  foreach(){
  	omega_out[] = omega[]*cs[];
  }

  view (fov = 21.0, tx = -0.5, ty = -0.5, width = 1200, height = 1200);
  clear();
  draw_vof ("f", lc = {0.0,0.0,0.0}, lw=1);
  draw_vof("cs", "fs",filled=-1);
  squares ("u.x", map = cool_warm, min = -0.1, max = 0.1);
  sprintf(timestring, "t=%2.02f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/HorizontalVelocity.mp4");

  view (fov = 21.0, tx = -0.5, ty = -0.5, width = 1200, height = 1200);
  clear();
  draw_vof ("f", lc = {0.0,0.0,0.0}, lw=1);
  draw_vof("cs", "fs",filled=-1);
  squares ("u.y", map = cool_warm, min = -0.1, max = 0.1);
  sprintf(timestring, "t=%2.02f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/VerticalVelocity.mp4");
  
  view (fov = 21.0, tx = -0.5, ty = -0.5, width = 1200, height = 1200);
  clear();
  draw_vof ("f", lc = {0.0,0.0,0.0}, lw=1);
  draw_vof("cs", "fs",filled=-1);
  squares ("f+cs", spread = -1, linear = true, map = cool_warm );
  sprintf(timestring, "t=%2.02f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/FluidPhases.mp4");

  view (fov = 21.0, tx = -0.5, ty = -0.5, width = 1200, height = 1200);
  clear();
  squares ("f+cs", spread = -1, linear = true, map = cool_warm );
  draw_vof ("f", lc = {0.0,0.0,0.0}, lw=1);
  draw_vof("cs", "fs",filled=-1);
  sprintf(timestring, "t=%2.02f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/FluidPhases_Background.mp4");
  
  view (fov = 21.0, tx = -0.5, ty = -0.5, width = 1200, height = 1200);
  clear();
  squares ("f+cs", spread = -1, linear = true, map = cool_warm );
  draw_vof ("f", lc = {0.0,0.0,0.0}, lw=1);
  draw_vof("cs", "fs",filled=-1);
  cells(lw = 1.0);
  sprintf(timestring, "t=%2.02f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/FluidPhasesGrid.mp4");

  view (fov = 21.0, tx = -0.5, ty = -0.5, width = 1200, height = 1200);
  clear();
  draw_vof ("f", lc = {0.0,0.0,0.0}, lw=1);
  draw_vof("cs", "fs",filled=-1);
  squares ("omega", map = cool_warm, min = -1.0, max = 1.0 );
  sprintf(timestring, "t=%2.02f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Vorticity.mp4");
  
  view (fov = 21.0, tx = -0.5, ty = -0.5, width = 1200, height = 1200);
  clear();
  draw_vof ("f", lc = {0.0,0.0,0.0}, lw=1);
  draw_vof("cs", "fs",filled=-1);
  squares ("omega_out", map = cool_warm, min = -1.0, max = 1.0 );
  sprintf(timestring, "t=%2.02f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Vorticity_Out.mp4");

}

event saveInterfaces (t = 0.0; t += 1.0; t <= tEnd) {

    char nameInterfaces1[200];

    sprintf(nameInterfaces1,"Interfaces/interfacesLiquid-%0.1f.dat",t);

    FILE * fp1 = fopen(nameInterfaces1, "w");
    output_facets (f, fp1);	
    fclose(fp1);
}

event loginterface (t += 0.1) {

    scalar posX[],posY[];
    position (f, posX, {1,0}); // (1,0) indicates the unit vector in the x-direction
    position (f, posY, {0,1}); // (0,1) indicates the unit vector in the y-direction

    fprintf(fp_interface, "%i %g %1.4f %1.4f %1.4f %1.4f %1.4f \n", i, t, statsf(f).sum, statsf(posX).min, statsf(posX).max, statsf(posY).min, statsf(posY).max);
    fflush(fp_interface);
}
  
event logvol (t += 0.1) {

    norm svol = normf(f);
    fprintf (fp_vol, "%d %g %g %g %g %g\n", i, t, svol.avg, svol.rms, svol.max, svol.volume);
    fflush(fp_vol);
}

event logstats (t += 0.1) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}
