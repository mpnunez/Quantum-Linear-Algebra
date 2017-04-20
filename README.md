# Quantum-Linear-Algebra
Numerical solution to the quantum harmonic oscillator (QHO) using Hartree-Fock and a plane wave basis set

### Description
Matlab code which uses a GUI. The user chooses the number of electrons in the system
and the size of the basis set (i.e. the number of plane waves). The user can choose
among several different graphs to display the solution. For example, the densities
computed numerically versus analytically are compared. Electron density is computed
using the many-electron wavefunction, as computed using a Slater determinant.

### History
Originally my final project for Spring 2013 CHEG827, University of Delaware, 
Department of Chemical and Biomolecular Engineering

Currently implmenting a GUI which will
- solve the equation with different parameters
- Display various parts of the solution
- Display the CPU time required to compute various parts
- Will handle not only 1 electron, but several with Hatree-Fock

See Report.pdf for details about the theory.