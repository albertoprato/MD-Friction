# Classical Molecular Dynamics & Friction Tensor Calculation

### Contents

1.  [Important Information](#1-important-information)
2.  [Features](#2-features)
3.  [Prerequisites](#3-prerequisites)
4.  [Program guide](#4-program-guide)

---

### 1. Important information

This Fortran program simulates the classical molecular dynamics of a solute consisting of four spheres immersed in a viscous solvent. The simulation uses Lennard-Jones potentials and integrates numerically the equations of motion using the Velocity Verlet algorithm.

The goal of this software is to compute the **friction tensor** of the solute as a function of time by analyzing the time-autocorrelation function of the forces experienced by the solute particles.

### 2. Features

* **Force Field**: Lennard-Jones potential with cutoff and Minimum Image Convention.
    * Solvent-Solvent interactions.
    * Solute-Solvent interactions.
* **Integration**: Velocity Verlet algorithm for stable time integration.
* **Minimization**: Conjugate Gradient (Polak-Ribiere) method to relax the initial solvent structure.
* **Initialization**: Maxwell velocity distribution
* **Analysis**: Computation of the friction tensor via Fast Fourier Transforms (FFTW3).

### 3. Prerequisites

To compile and run this simulation, you need:

1.  **Make**: GNU Make build tool.
2.  **GFortran**: The GNU Fortran compiler.
3.  **FFTW3**: The "Fastest Fourier Transform in the West" library.

### 4. Program guide

To understand how to use this Fortran program, please read the User_Guide.pdf manual located in the doc folder.
