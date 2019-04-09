%-------------------------------------------------------------------------%
%
%           Scattering of by a human head 
%
%    The permittivities within a head are extremely high and have a large
%    range - between 0 and 100. Therefore the preconditioning challenge is
%    somewhat different to the high-frequency low-contrast ice crystal 
%    setting.
%    This script aims to test the effectiveness of the 2-level circulant 
%    preconditioner for this challenging scenario. We must use some
%    permittivity averaging in order to set up the preconditioner.
%    
%    S.P.Groth 4th Apr 2019
%
%-------------------------------------------------------------------------%
close all; 
clear;
format compact;
%% Add the directory path
addpath(genpath('piecewise_constant'))

%% Geometrical parameters
% load('meshes/Head_model_Billie.mat')
load('meshes/Head_model_Duke.mat')
% load('meshes/Body_model_Billie.mat')
r = RHBM.r;
epsilon_r = RHBM.epsilon_r;
sigma_e = RHBM.sigma_e;
idx = RHBM.idxS;
lambda_ext = 550e-6; 

epsilon_r(idx) = epsilon_r(idx);



form=2;
freq = 300*1e6;
EMconstants;

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
% keyboard
%% Define incident field (y-polarized plane wave travelling in x-direction)
Eo = [0, 0, 1];     % z-polarization
dInc = [0, 1, 0];   % direction
% ko = 2*pi/lambda_ext;
kvec = ko * dInc;
[Einc, ~] = PlaneWave_Excitation(r,kvec,omega_mu,Eo);

%% Generate operator
[fN,opToeplitz] = getOPERATORS_SAM(r,freq,'N',[],'DEMCEM');

%% Establish DDA matrix-vector product function
Mr = epsilon_r;
Mc = epsilon_r - 1.0;

% get the positions of the non-air voxels in the 3D grid
idxS = find(abs(Mc(:)) > 1e-12); % these are the indexes of the non-air voxel positions
idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions for 3 Cartesian components

if form==2
    Mc(:) = Mc(:)./Mr(:);
    % notice that Vrhs is only defined for non-air components
    Mr(:) = Mr(:)./Mr(:);
    
end

% get the voxel side
dx = r(2,1,1,1) - r(1,1,1,1);

% Gram matrix (value)
Gram = dx^3; % volume of the voxel

% get the positions of the non-air voxels in the 3D grid
idx3 = [idx; nD+idx; 2*nD+idx]; % the vector of non-air positions for 3 Cartesian components

% set the multiplier to convert E fields to the solver rhs (currents)
tau = 1j*omega*eo*Mc;
tau3 = [tau(:); tau(:); tau(:)]; % 3 Cartesian components in vector form

% get the RHS for the solver
% notice that only the non-air positions in Cartesian components are used
Vrhs = Gram.*tau3(idxS3).*Einc(idxS3);

if form==2
    if L==1
        Vrhs(:) = Vrhs(:)./[squeeze(Mr(idxS)); squeeze(Mr(idxS)); squeeze(Mr(idxS))];
    else
        Vrhs(:) = Vrhs(:)./[Mr(idxS); Mr(idxS); Mr(idxS)]; % notice that Vrhs is only defined for non-air components
    end
end

fACPU   = @(J)mv_AN_const(J, fN, Mr, Mc, Gram, 'notransp', idxS3, 0);

tol = 1e-12;
 % Solve without preconditioner
 tini=tic;
%  [x,flag,relres,iter,resvec] = gmres(@(J)fACPU(J), Vrhs,500, tol, 1);
[x,flag,relres,iter,resvec] = bicgstab(@(J)fACPU(J), Vrhs, tol, 2000);
 tend=toc(tini);
 nIts_0 = length(resvec);
 fprintf('No preconditioner. Solve time = %.2f [sec] \n',tend)
 fprintf('Iteration count = %d \n',nIts_0);
 
 Je = zeros(L,M,N,3);
Je(idxS3) = x ;

DOF=3*L*M*N

[Eout] = E_field_Nop_const(Je,fN,Gram,freq,Einc);

figure
pcolor(real(rot90(Eout(:,:,round(N/2),3))))
shading interp
colormap(bluewhitered)
axis image


%% Preconditioner assembly
% Compute 1-level circulant approximation to Toeplitz operator
tic
[circ_N,circ_L_opToep] = circ_1_level(opToeplitz, L, M, N, 0, 'off');

% Compute 2-level circulant approximation to Toeplitz operator
[circ_M_opToep,circ_2_N] = circ_2_level(circ_L_opToep,L,M,N);
disp('Circulant approximation DDA Toeplitz operator');   
toc

% % Average permittivities along x direction
% hey = squeeze(sum(Mc,1)./L);
% Mc_av_x = ones(L,M,N);
% for i=1:L
%     Mc_av_x(i,:,:) = hey;
% end
% % Now average along the y direction
% kk = squeeze(sum(Mc_av_x,2)./M);
% Mc_av_xy = ones(L,M,N);
% for i=1:M
%     Mc_av_xy(:,i,:) = kk;
% end
% face_2=squeeze(Mc_av_xy(1,1,:));
% face3_2 = [face_2(:);face_2(:);face_2(:)];
% face3_2mat = diag(face3_2);

% Average permittivities along x direction
hey = squeeze(sum(epsilon_r,1)./L);
er_av_x = ones(L,M,N);
for i=1:L
    er_av_x(i,:,:) = hey;
end

temp0 = squeeze(er_av_x(1,:,:));
face=(temp0-1)./temp0;
face3 = [face(:);face(:);face(:)];
face3mat = diag(face3);

% Now average along the y direction
kk = squeeze(sum(er_av_x,2)./M);
er_av_xy = ones(L,M,N);
for i=1:M
    er_av_xy(:,i,:) = kk;
end
temp=squeeze(er_av_xy(1,1,:));
face_2 = (temp-1)./temp;
face3_2 = [face_2(:);face_2(:);face_2(:)];
face3_2mat = diag(face3_2);

tic;
for i=1:L
    for j=1:M
        circ_2_temp = Gram*eye(3*N)-face3_2mat*circ_2_N{i,j};
        circ_2_inv{i,j} = inv(circ_2_temp);
    end
end
disp('Inversion of 2-level circulant preconditioner');
toc

%  tic
% circ_inv=cell(L);
% parfor i=1:L
%     circ_temp = Gram*eye(3*M*N)-face3mat*circ_N{i};
%     circ_inv{i} = inv(circ_temp);
% end
% disp('Inversion of 1-level circulant preconditioner');
% toc

% MVP function for precondioner inverse
prec = @(J) mvp_circ_2_level(circ_2_inv,J,L,M,N,idx3);

% prec1 = @(J) chan_mvp_idx(circ_inv,J,L,M,N,idx3);

tini1_gmres=tic;
% [vsol1_gmres,~,~,~,resvec1_gmres] = gmres(@(J)fACPU(J), Vrhs, 2000,...
%     tol, 1, @(J)prec(J));
[vsol1_gmres,~,~,~,resvec1_gmres] = bicgstab(@(J)fACPU(J), Vrhs, tol, 2000,...
    @(J)prec(J));
tend1_gmres=toc(tini1_gmres);
nIts_1_gmres = length(resvec1_gmres);
fprintf('GMRES WITH 2-level preconditioner. Solve time = %.2f [sec] \n',tend1_gmres)
fprintf('Iteration count = %d \n',nIts_1_gmres);

% tini1_gmres=tic;
% [vsol1_gmres,~,~,~,resvec1_gmres] = gmres(@(J)fACPU(J), Vrhs, 2000,...
%     tol, 1, @(J)prec1(J));
% tend1_gmres=toc(tini1_gmres);
% nIts_1_gmres = length(resvec1_gmres);
% fprintf('GMRES WITH 1-level preconditioner. Solve time = %.2f [sec] \n',tend1_gmres)
% fprintf('Iteration count = %d \n',nIts_1_gmres);

