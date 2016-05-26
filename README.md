# Quantum-Linear-Algebra
Solution of the Quantum harmonic oscillator using a plane wave basis set

Final project for Spring 2013 CHEG827

Shows how to formulate a basic quantum chemistry problem as an eigenvalue problem. 
Computational time is saved by using eigs instead of eig in Matlab. Eigs uses Krylov
methods which only seek the important eigenvalues we need.

See Report.pdf for details.