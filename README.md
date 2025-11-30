# 2D Liquid Bridge Dynamics

Direct numerical simulation code infrastructure for the two-dimensional computational model of a liquid bridge created between two cylinders held at a fixed distance between them, and subsequently allowed to fall under gravity. The key challenge is including contact angle dynamics over embedded solid geometries, as only recently (2024) [implemented in Basilisk](http://basilisk.fr/sandbox/tavares/) and expanded upon in [the associated publication](https://www.sciencedirect.com/science/article/pii/S004579302400149X). The code supports the active collaboration on this topic with the Department of Mathematical Sciences and Department of Chemical Engineering at Loughborough University.

<img width="100%" alt="LiquidBridge" src="https://github.com/user-attachments/assets/77e9782f-d896-41fb-a2d3-65d715d8a5cf" />


---

## 📌 Highlights

- 2D Basilisk-based solver for liquid bridge dynamics between cylinders  
- Includes contact angle dynamics on embedded solids (contact-embed implementation)  
- Useful for parameter sweeps in geometry, spacing and fluid properties  
- Lightweight driver scripts to run example cases and organize outputs

---

## 🛠️ Installation

### 1. Requirements

- [Basilisk](http://basilisk.fr/) (compiled with `qcc`)  
- C compiler (gcc/clang)  
- Visualization tools for post-processing (e.g., ffmpeg, imagemagick, gnuplot) — recommended:
  ```bash
  sudo apt install ffmpeg imagemagick gnuplot
  ```

### 2. Clone the repository

```bash
git clone https://github.com/rcsc-group/LiquidBridge2D
cd LiquidBridge2D
```

### 3. Install Basilisk

Follow Basilisk’s installation instructions: http://basilisk.fr/src/INSTALL

### 4. Running the example case

A convenient shell script is provided to run example cases and collect outputs:

```bash
sh DriverCode/run_bridge_Fr.sh
```

The script invokes the driver code and places outputs into organized folders for analysis and visualization.

---

## ⚙️ Key Simulation Options

Edit parameters in the driver/source to control:

- Effective nondimensional numbers (e.g., Froude or other scalings) via the driver code  
- Cylinder spacing and radii (in the driver/source)  
- Resolution / grid refinement (in the driver/source)  
- Output frequency and visualization options (inside the driver code)

Key implementation notes:
- DriverCode/contact-embed.h contains the embedded-solid / contact-angle routines used by the solver.
- The main driver is DriverCode/liquid_bridge_Fr.c.
- run_bridge_Fr.sh wraps execution and output organization.

---

## 📁 Folder Structure

```bash
.
├── DriverCode/                  # Main driver scripts and code for simulations
│   ├── contact-embed.h          # Embedded-solid & contact-angle implementation
│   ├── liquid_bridge_Fr.c       # Main simulation driver/source (Basilisk)
│   ├── run_bridge_Fr.sh         # Shell script to run example cases and manage outputs
│
├── LICENSE                      # License information
└── README.md                    # Project documentation
```

---

## 📊 Outputs

Typical outputs produced by the driver/script:

- Interface snapshots and data files for post-processing  
- Time series / summary files with diagnostic quantities
- Optional animations (if visualization routines are enabled)

Output layout and file naming are managed by the run script; inspect DriverCode/run_bridge_Fr.sh for details.

---

## 🧑 Contributing

Contributions welcome:

- Fork the repo  
- Open issues to report bugs or request features  
- Submit pull requests with improvements, tests, or documentation updates
