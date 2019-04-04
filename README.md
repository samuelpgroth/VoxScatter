# VoxScatter

Matlab repository for computing electromagnetic scattering by dielectric particles.

The code solves the electric volume integral equation via a choice of two (similar) techniques:

1. The Discrete Dipole Approximation (DDA) based on: 
B T Draine and P J Flatau. Discrete-dipole approximation for scattering calculations. JOSA A, 11(4):1491–1499, 1994.

2. The Galerkin Method of Moments (MoM) based on:
A G Polimeridis, J Fernandez Villena, L Daniel, and J K White. Stable FFT-JVIE solvers for fast analysis of highly
inhomogeneous dielectric objects. Journal of Computational Physics, 269:280–296, 2014. This uses code from the
repository https://github.com/thanospol/MARIE

Both approaches use a voxelized (uniform) discretization of the particle. This enables the acceleration of 
matrix-vector products with the fast-Fourier transform (FFT). The second approach has better conditioning properties, especially for large refractive indices.

## Circulant preconditioning

The convergence of the iterative solves are hugely accelerated by the use of preconditioners based on the circulant 
approximation of the system matrix. Efficient implementations of the 1- and 2-level circulant preconditioners of 
Chan and Olkin are included. See our preprint https://arxiv.org/abs/1903.09802 for details.
