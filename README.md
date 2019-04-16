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
Chan and Olkin are included. See our preprint https://arxiv.org/abs/1903.09802 for details. The results in this preprint can be generated from the script main_DDA.m

## Current example scripts
main_DDA.m performs scattering by hexagonal plates. We demonstrate the performance of circulant preconditioning by solving the linear system four times: twice without preconditioning (via BiCG-Stab, GMRES) and twice with preconditiong (via BiCG-Stab, GMRES). Should agree with results presented in our preprint https://arxiv.org/abs/1903.09802

main_koch.m is the same as above but for Koch fractal snowflakes.

main_koch_prec_only.m is the same as main_koch.m but with only one solve - with preconditioning via BiCG-Stab (the fastest). This produces pretty pictures, for example:
![koch_scatter](https://user-images.githubusercontent.com/13260045/56199494-efd83280-6034-11e9-8966-5851276a31f5.png)

## To do:
1. Add the volume integral equation (rigorous implementation) as an option. Show comparision with "standard" DDA.

2. Far-field scattering routine.

3. Random orientation.
