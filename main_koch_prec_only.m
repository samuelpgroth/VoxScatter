%-------------------------------------------------------------------------%
%
%           Scattering by a single particle with DDA 
%
%    Compute the field scattered by a plain wave incident upon a dielectric
%    obstacle using the Discrete Dipole Approximation (DDA).
%
%    The linear system is solved iteratively with BiCGStab and GMRES, first
%    without any preconditioning, then with a 2-level circulant
%    preconditioner.
%
%    Compare the results of the current setup to Table 3 in the paper
%    draft.
%    
%    S.P.Groth 6th Feb 2019
%
%-------------------------------------------------------------------------%
close all; 
clear;
format compact;
addpath(genpath('src'))

%% Geometrical parameters
sizeParam = 20;  % size parameter
nPerLam = 10;    % number of voxels per interior wavelength

% Refractive index of scatterer (real and imaginary parts)
refRe = 1.5;
refIm = 0;

%% Generate voxel coordinates and vizualize
refInd = refRe + 1i*refIm;
% [r, idx, res, P, lambda_ext, lambda_int] = ...
%     geometry(geom, refInd, sizeParam, nPerLam);

% Need to neaten this Koch example and place into geometry.m
lambda_ext = 2*pi/sizeParam; % exterior wavelength of incident wave
lambda_int = lambda_ext/real(refInd); % interior wavelength
aspectRatio = 1/20;
koch_snowflake;
idx = find(in);

[L,M,N,~]=size(r);
nD = L * M * N;   % number of voxels in simulation domain

% plot geometry
xd = r(:,:,:,1);
yd = r(:,:,:,2);
zd = r(:,:,:,3);

figure(1);
plot3(xd(idx), yd(idx), zd(idx), 's');
axis equal;
grid on;

%% Define incident field (y-polarized plane wave travelling in x-direction)
Eo = [0, 0, 1];     % z-polarization
dInc = [1, 0, 0];   % direction
ko = 2*pi/lambda_ext;
kvec = ko * dInc;
[Einc] = PlaneWavePointWise(Eo, kvec, r);

%% Generate operator
% Perform nearby quadrature more accurately 
nearby_quad = 'off'; % 'off' = standard DDA, 'on' ~= collocation VIE

% Generate Toeplitz DDA operator and its circulant embedding
[operator_circ, operator_toep, alpha_LDR] = ...
    getOPERATOR_DDA(r, ko, refInd, kvec, Eo, nearby_quad);

%% Establish DDA matrix-vector product function
id = ones(L, M, N);   % identity
chi = zeros(L, M, N);
chi(idx) = 1;       
dx = r(2,1,1,1) - r(1,1,1,1);
dV = dx ^ 3;   % voxel volume

% get the positions of the non-air voxels in the 3D grid
idx3 = [idx; nD+idx; 2*nD+idx]; % the vector of non-air positions for 3 Cartesian components

% Matrix-vector product
MVP_DDA = @(J)MVP(J, operator_circ, id, dV*chi, 1/alpha_LDR, 'notransp', idx3, 0);

%% Preconditioner assembly
% Compute 1-level circulant approximation to Toeplitz operator
tic
[circ_N,circ_L_opToep] = circ_1_level(operator_toep, L, M, N, 0, 'off');

% Compute 2-level circulant approximation to Toeplitz operator
[circ_M_opToep,circ_2_N] = circ_2_level(circ_L_opToep,L,M,N);
disp('Circulant approximation DDA Toeplitz operator');   
toc

% Invert diagonal blocks of circulant approximation
tic
parfor i=1:L
    for j=1:M
        circ_2_inv{i,j} = inv(1/alpha_LDR*eye(3*N)-dV*circ_2_N{i,j});
    end
end
disp('Inversion of circulant preconditioner')
toc

% MVP function for precondioner inverse
prec = @(J) mvp_circ_2_level(circ_2_inv,J,L,M,N,idx3);

%% Solve system with PRECONDITIONED iterative method
tol = 1e-5;
Vrhs = Einc(idx3);
tini1=tic;
[vsol1,~,~,~,resvec1] = bicgstab(@(J)MVP_DDA(J), Vrhs, tol, 2000,...
    @(J)prec(J));
tend1=toc(tini1);
nIts_1 = length(resvec1);
fprintf('BICG WITH preconditioner. Solve time = %.2f [sec] \n',tend1)
fprintf('Iteration count = %d \n',nIts_1);

%% Gather and visualize total (incident + scattered) field
Pout = zeros(L,M,N,3);  % Polarizations
Eout = zeros(L,M,N,3);  % Electric field
Pout(idx3) = vsol1;
CHI = (refInd^2 - 1) / (4 * pi);
Eout(idx3) = vsol1./CHI;

% Evaluate field throughout computation domain
idxAll = (1:3*nD)';
MVP_eval = @(J)MVP(J, operator_circ, id, dV*id, 1/alpha_LDR, ...
    'notransp', idxAll, 0);
temp = MVP_eval(Pout(:));

Etemp = zeros(L,M,N,3);
Etemp(idxAll) = temp;
E = Einc - Etemp + Eout;
figure
pcolor(xd(:,1,2),yd(1,:,2),rot90(real(E(:,:,3,3))))
shading interp
axis image
colormap bluewhitered
colorbar
hold on
plot(P(:,1),P(:,2),'k','LineWidth',2);
set(gcf,'color','w')