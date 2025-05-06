// From http://basilisk.fr/sandbox/tavares/

(const) scalar contact_angle;

static inline coord normal_contact (coord ns, coord nf, double angle)
{
  assert (dimension == 2); // fixme: 2D only for the moment
  coord n;
  if (- ns.x*nf.y + ns.y*nf.x > 0) { // 2D cross product
    n.x = - ns.x*cos(angle) + ns.y*sin(angle);
    n.y = - ns.x*sin(angle) - ns.y*cos(angle);
  }
  else {
    n.x = - ns.x*cos(angle) - ns.y*sin(angle);
    n.y =   ns.x*sin(angle) - ns.y*cos(angle);
  }
  return n;
}

void reconstruction_contact (scalar f, vector n, scalar alpha)
{

  reconstruction (f, n, alpha);
  
    foreach()
    if (cs[] < 1. && cs[] > 0. && f[] < 1.  && f[] > 0.) {
      coord ns = facet_normal (point, cs, fs);
      normalize (&ns);
      coord nf;
      foreach_dimension()
	nf.x = n.x[];
      coord nc = normal_contact (ns, nf, contact_angle[]);
      foreach_dimension()
	n.x[] = nc.x;
      alpha[] = line_alpha (f[], nc);
    }
  boundary ({n, alpha});  
}

event contact (i++)
{
  vector n[];
  scalar alpha[];
  reconstruction_contact (f, n, alpha);
  foreach() {
    if (cs[] == 0.) {
      double fc = 0., sfc = 0.;
      coord o = {x, y, z};
      foreach_neighbor()
	if (cs[] < 1. && cs[] > 0. && f[] < 1.  && f[] > 0.) {
	  double coef = cs[]*(1. - cs[])*f[]*(1. - f[]);
	  sfc += coef;
	  coord nf;
	  foreach_dimension()
	    nf.x = n.x[];
	  coord a = {x, y, z}, b;
	  foreach_dimension()
	    a.x = (o.x - a.x)/Delta - 0.5, b.x = a.x + 1.;
	  fc += coef*rectangle_fraction (nf, alpha[], a, b);
	}
      if (sfc > 0.)
	f[] = fc/sfc;
    }
  }
  boundary ({f});
}
