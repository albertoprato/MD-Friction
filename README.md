## Classical Molecular Dynamics & Friction Tensor Calculation

**Author**:
  * Alberto Prato - [alberto.prato@studenti.unipd.it](mailto:alberto.prato@studenti.unipd.it)

---

### Contents

1.  [Information](#1-information)
2. [Repository Structure](#2-repository-structure)
3.  [Prerequisites](#3-prerequisites)
4.  [Program guide](#4-program-guide)

---

### 1. Information

This Fortran program simulates the classical molecular dynamics of a solute consisting of four spheres immersed in a viscous solvent. The simulation uses Lennard-Jones potentials and numerically integrates the equations of motion using the Velocity Verlet algorithm.

The goal is to compute the **friction tensor** of the solute as a function of time. This is achieved by calculating the time-autocorrelation function of the forces experienced by the solute particles during the simulation.

For further information regarding the physical context and underlying equations, please refer to the `Methodology_and_Implementation.pdf` located in the `doc` folder.

---

### 2. Repository Structure

The repository is organized into three main directories:

* **`verlet-LJ/`**: contains the Fortran source code for the molecular dynamics simulation, including the Makefile for compilation.
* **`python_tool/`**: Includes a jupyter notebook designed to visualize the output data.
* **`doc/`**: Documentation folder containing:
    * **`User_Guide.pdf`**: a practical manual explaining how to compile, configure, and run the software.
    * **`Scientific_Report.pdf`**: a detailed report covering the theoretical chemical-physical background and coding strategies.

---

### 2. Prerequisites

To compile and run this simulation, the following tools and libraries are required: 

1.  **Make**: The GNU Make build tool.
2.  **GFortran**: The GNU Fortran compiler.
3.  **FFTW3**: The "Fastest Fourier Transform in the West" library.

---

### 3. Program guide

To understand how to use this Fortran program, please read the `User_Guide.pdf` manual located in the `doc` folder.
