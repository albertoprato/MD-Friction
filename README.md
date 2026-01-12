## Classical Molecular Dynamics & Friction Tensor Calculation

* Alberto Prato - alberto.prato@studenti.unipd.it

---

### Contents

1.  [Information](#1-information)
2.  [Prerequisites](#2-prerequisites)
3.  [Program guide](#3-program-guide)

---

### 1. Information

This Fortran program simulates the classical molecular dynamics of a solute consisting of four spheres immersed in a viscous solvent. The simulation uses Lennard-Jones potentials and integrates numerically the equations of motion using the Velocity Verlet algorithm.

The goal is to compute the **friction tensor** of the solute as a function of time through the time-autocorrelation function of the forces experienced by the solute particles.

---

### 2. Prerequisites

To compile and run this simulation, you need:

1.  **Make**: GNU Make build tool.
2.  **GFortran**: The GNU Fortran compiler.
3.  **FFTW3**: The "Fastest Fourier Transform in the West" library.

---

### 3. Program guide

To understand how to use this Fortran program, please read the User_Guide.pdf manual located in the doc folder.
