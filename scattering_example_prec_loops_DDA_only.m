%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Scattering by a single particle
%
%  Calculate the scattering matrix entries for an unpolarized plane wave
%  travelling in the x-direction incident upon a particle in fixed
%  orientation.
%
%  Very much a development version - email me at samgroth@mit.edu with
%  questions/comments/suggestions.
%
%  The only portions you need to worry about are the "Parameter entry" at
%  the top and the "Outputs" at the bottom.
%
%  S.P.Groth 13/2/17
%
%-------------------------------------------------------------------------%
close all
clear;
format compact
% format long
%% Add the corresponding path
addpath(genpath('piecewise_constant'))

% Add a test comment

SIZE_PARAM = [10 20 30 40 60 80 100];
REF_RE = [1.2 1.6 1.8 2];
REF_NAME = {'1p2','1p6','1p8','2'};

for i_param = 6:length(SIZE_PARAM)
    for i_ref = 1:1%length(REF_RE)

%% Parameter entry

sizeParam = SIZE_PARAM(i_param);%3.788399944423817;  % size paramter ko*a
nPerLam = 10;    % number of voxels per wavelength (10 for rough solution, 20 for one percent or better)
geom = 'hex'; % current choice: hex, sphere, spheroid (please ask me to add more, if required)

% Refractive of ice
refRe = REF_RE(i_ref);      % real part
refIm = 0;%2.289e-9;   % imaginary part

% Gometry parameters (alter appropriate geometry)
if strcmp(geom,'spheroid')==1
    a = 1;   % a is radius in x- and y-directions
    c = 0.5; % c is radius in z-direction. c<a => oblate spheroid
    dom_x = 2*a;
    dom_y = 2*a;
    dom_z = 2*c;
    P = pi*a*c;  % projected  cross-section
elseif strcmp(geom,'sphere')==1
    a = 1;
    c = 1;
    dom_x = 2*a;
    dom_y = 2*a;
    dom_z = 2*a;
    P = pi*a^2;  % projected cross-section
elseif strcmp(geom,'hex')==1
    a = 1;
    b = sqrt(3)/2*a; % half-height of the hexagonal face
    aspectRatio = 1/5; % ratio of column's height to it's width
    dom_x = 2*a;
    dom_y = 2*b;
    dom_z = a*aspectRatio;
    P = dom_y*dom_z;
elseif strcmp(geom,'cube')==1
    a = 1;
    dom_x = 2*a;
    dom_y = 2*a;
    dom_z = 2*a;
    P = 4*a^2;
end

refInd = refRe+1i*refIm;

%% Define some standard EM quantities
lambda_ext = 2*pi*a/sizeParam; % exterior wavelength of incident wave
co = 299792458; % speed of light
freq = co/lambda_ext;
EMconstants;
% This code take the permittivity, not refractive index, so we square the
% above and take the real and imaginary parts.
e_r = real(refInd^2);  % real part
s_e = imag(refInd^2);  % imaginary part

%% DISCRETIZATION (into voxels)
% Preferential dimension:
%this is the direction in which we ensure that the voxels match the
%geometry
h_pref = dom_x; % here it's the height (z-direction)
lambda_int = lambda_ext/sqrt(e_r); % interior wavelength
% figure out resolution
res_temp = lambda_int/nPerLam; % provisional resolution
N = ceil(h_pref/res_temp);
res = h_pref/N;

% Generate voxels
[r] = generatedomain_new(res,dom_x/2,dom_y/2,dom_z/2);
[L,M,N,~]=size(r);

% Find the indices corresponding to those voxels inside the shape
if strcmp(geom,'hex')==1    % Hexagonal column
    % Each line defines a sloping face of the hexagon (the top and bottom
    % ones are easier)
    l1 = @(r) (r(:,:,:,2)-sqrt(3)*(a-r(:,:,:,1)));
    l2 = @(r) (r(:,:,:,2)-sqrt(3)*(r(:,:,:,1)+a));
    l3 = @(r) (r(:,:,:,2)-sqrt(3)*(-a-r(:,:,:,1)));
    l4 = @(r) (r(:,:,:,2)-sqrt(3)*(r(:,:,:,1)-a));

    point = (l1(r)<0).*(l2(r)<0).*(l3(r)>0).*(l4(r)>0);
    idx = find(point);  % indices

elseif strcmp(geom,'spheroid')==1 || strcmp(geom,'sphere')==1
    % Spheroid with radius a in x- and y- directions, and radius c in
    % z-direction
    spheroidFun = @(r) (r(:,:,:,1).^2+r(:,:,:,2).^2)./a^2 + (r(:,:,:,3).^2) ./ c^2;
    idx = find(spheroidFun(r)<=1);
elseif strcmp(geom,'cube')==1
    idx = 1:L*M*N;
end

% plot geometry
xd = r(:,:,:,1);
yd = r(:,:,:,2);
zd = r(:,:,:,3);

numVox = L*M*N;
numVoxIce = length(idx);

figure(1);
plot3(xd(idx), yd(idx), zd(idx), 's');
axis equal;
grid on;

epsilon_r = ones(L,M,N);
epsilon_r(idx) = e_r;

sigma_e = zeros(L,M,N);
sigma_e(idx) = s_e;

epsilon_r = epsilon_r - 1j*sigma_e;  % note minus sign, following the engineering literature we use -1j instead of +1i

% get the voxel side
dx = r(2,1,1,1) - r(1,1,1,1);

% Gram matrix (value)
Gram = dx^3; % volume of the voxel

%% Sub-voxel smoothing

av=0;

if av==1

    n_f=5;

    res_f = res/n_f;

    [r_f] = generatedomain_new(res_f,dom_x/2,dom_y/2,dom_z/2);
    [L_f,M_f,N_f,~]=size(r_f);

    xd_f = r_f(:,:,:,1);
    yd_f = r_f(:,:,:,2);
    zd_f = r_f(:,:,:,3);

    idx_f = find( xd_f.^2+yd_f.^2+zd_f.^2<=a^2);

    epsilon_r_f = ones(L_f,M_f,N_f);

    epsilon_r_f(idx_f) = e_r;


    for ix=1:L
        for iy=1:M
            for iz=1:N
                temp_er=epsilon_r_f((ix-1)*n_f+1:ix*n_f,(iy-1)*n_f+1:iy*n_f,...
                    (iz-1)*n_f+1:iz*n_f);

                epsilon_r(ix,iy,iz) = sum(sum(sum(temp_er)))./(n_f^3);
            end
        end
    end
end

% keyboard

% -------------------------------------------------------------------------
%                  Generate circulants
% -------------------------------------------------------------------------

% compute the circulants
% [fN,opToep] = getOPERATORS_SAM(r,freq,'N',[],'DEMCEM');

%% First polarization solve (z-polarized plane wave)

% Define excitation - a plane wave
Eo = [0,0,1]; % z-polarized incident E-field
k = ko * [1,0,0]; % travelling in +ve x-direction
[Einc, ~] = PlaneWave_Excitation(r,k,omega_mu,Eo);

% Solve system
%% Assemble the entire matrix to look at e-values and condition number
% Compute the relative permittivity and suceptibility
Mr = epsilon_r;
Mc = epsilon_r - 1.0;

% Domain dimensions
[L,M,N,~] = size(r);
nD = L*M*N; % number of variables in the system

% get the positions of the non-air voxels in the 3D grid
idxS = find(abs(Mc(:)) > 1e-12); % these are the indexes of the non-air voxel positions
idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions for 3 Cartesian components

% form = 2; % Always solve JVIE 2.
% if form==2
%     Mc(:) = Mc(:)./Mr(:);
%     % notice that Vrhs is only defined for non-air components
%     Mr(:) = Mr(:)./Mr(:);
%
% end
%
% tau = 1j*omega*eo*Mc;
% tau3 = [tau(:); tau(:); tau(:)]; % 3 Cartesian components in vector form

% get the RHS for the solver
% notice that only the non-air positions in Cartesian components are used
% Vrhs = Gram.*tau3(idxS3).*Einc(idxS3);
%
% if form==2
%     if L==1
%         Vrhs(:) = Vrhs(:)./[squeeze(Mr(idxS)); squeeze(Mr(idxS)); squeeze(Mr(idxS))];
%     else
%         Vrhs(:) = Vrhs(:)./[Mr(idxS); Mr(idxS); Mr(idxS)]; % notice that Vrhs is only defined for non-air components
%     end
% end
%
% fACPU   = @(J)mv_AN_const(J, fN, Mr, Mc, Gram, 'notransp', idxS3, 0);
%
% tini2 = tic;
%
% tic
% [circ_N,circ_L_opToep] = level_1_parallel_func_N(opToep,L,...
%     M,N,Gram,0,'off');
% disp('Circulant approximation of N operator');
% toc


% [circ_M_opToep,circ_2_N] = circ_2_level_fast_test(circ_L_opToep,L,M,N);
% %
% for i=1:L
%     for j=1:M
%         circ_2_inv_T{i,j} = inv(Gram*eye(3*N)-(e_r-1)/e_r*circ_2_N{i,j});
%     end
% end
% disp('2-level assembly')
% t_2level = toc(tini2)
%
% prec2 = @(J) chan_2_mvp_for_parallel_idx(circ_2_inv_T,J,L,M,N,idxS3);
%
% tic
% [circ_2_inv_T] = level_2_parallel_func(opToep,Mc,Mr,L,M,N,Gram);
% disp('Second level assembly');
% toc


%--------------------------- Solve system --------------------------------%
% tol = 1e-5;
% % Solve without preconditioner
% tiniV0=tic;
% [vsolV0,~,~,~,resvecV0] = pgmres(@(J)fACPU(J), Vrhs,2000, tol, 1);
% %[vsolV0,~,~,~,resvecV0] = bicgstab(@(J)fACPU(J), Vrhs, tol, 2000);
% tendV0=toc(tiniV0);
% nIts_0 = length(resvecV0);
% fprintf('No preconditioner. Solve time = %.2f [sec] \n',tendV0)
% fprintf('Iteration count = %d \n',nIts_0);
%
% % Solve without preconditioner
% tiniV1=tic;
% [vsolV1,~,~,~,resvecV1] = pgmres(@(J)fACPU(J), Vrhs,2000, tol, 1,@(J)prec2(J));
% %[vsolV1,~,~,~,resvecV1] = bicgstab(@(J)fACPU(J), Vrhs, tol, 2000,@(J)prec2(J));
% tendV1=toc(tiniV1);
% nIts_0 = length(resvecV1);
% fprintf('Preconditioner. Solve time = %.2f [sec] \n',tendV1)
% fprintf('Iteration count = %d \n',nIts_0);

% keyboard

kvec=k;

%% Self-interaction Classius-Mossotti stuff

b1 = -1.8915316;
b2 = 0.1648469;
b3 = -1.7700004;
msqr = refInd.^2;
dcube = Gram;
d=dx;

% if nargin > 3  % we have polarization info
  a_hat = kvec/norm(kvec);
  e_hat = Eo/norm(Eo);
  S = 0;
  for j = 1:3
    S = S + (a_hat(j)*e_hat(j))^2;
  end
% else           % use randomly-oriented value; also for non-plane wave
%   S = .2;
% end

alpha_CM = 3*dcube/(4*pi)*(msqr - 1)./(msqr + 2); % Clausius-Mossotti
alpha_LDR = alpha_CM./(1 + (alpha_CM/dcube).*((b1+msqr*b2+msqr*b3*S)*(ko*d)^2-2/3*1i*ko^3*dcube));

%%

Mr = ones(L,M,N);
Mc = zeros(L,M,N);
Mc(idx) = 1;

nD = L*M*N;
idxS = find(abs(Mc(:)) > 1e-12); % these are the indexes of the non-air voxel positions
idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions for 3 Cartesian components

close all

%%
I = eye(3);

Toep = zeros(L,M,N,3);
R0 = squeeze(r(1,1,1,:));

nQuad=1;
[wG,xG] = gauss_1d(nQuad);
[XG,YG,ZG] = meshgrid(xG);
[XW,YW,ZW] = meshgrid(wG*0.5);

tIni = tic;

for i=1:L
    for j=1:M
        for k=1:N
            R1 = squeeze(r(i,j,k,:));
            rk_to_rj = R1-R0;
            rjk = norm(rk_to_rj); %sqrt(sum((r(jj,:)-r(kk,:)).^2))

%             if rjk<3*dx && rjk>1e-15
% %                 keyboard
%                 x_grid = R1(1) + dx/2 * XG;
%                 y_grid = R1(2) + dx/2 * YG;
%                 z_grid = R1(3) + dx/2 * ZG;
%
%                 temp=zeros(3,3);
%                 for iQ = 1:nQuad
%                     for jQ = 1:nQuad
%                         for kQ = 1:nQuad
%                             RQ = [x_grid(iQ,jQ,kQ);y_grid(iQ,jQ,kQ);...
%                                 z_grid(iQ,jQ,kQ)];
%
%
%                             rk_to_rj = RQ-R0;
%
%
%
%                             rjk = norm(rk_to_rj); %sqrt(sum((r(jj,:)-r(kk,:)).^2))
%                             rjk_hat = (rk_to_rj)/rjk;
%                             rjkrjk = rjk_hat*rjk_hat';
%                             %             keyboard
%
%                             Ajk = exp(1i*ko*rjk)/rjk*(ko^2*(rjkrjk - I) + (1i*ko*rjk-1)/rjk^2*(3*rjkrjk - I)); %Draine & Flatau
%                             temp = temp + Ajk.*XW(iQ,jQ,kQ).*...
%                                 YW(iQ,jQ,kQ).*ZW(iQ,jQ,kQ);
%
%                         end
%                     end
%                 end
%                 Toep(i,j,k,1) = temp(1,1);
%                 Toep(i,j,k,2) = temp(1,2);
%                 Toep(i,j,k,3) = temp(1,3);
%                 Toep(i,j,k,4) = temp(2,2);
%                 Toep(i,j,k,5) = temp(2,3);
%                 Toep(i,j,k,6) = temp(3,3);
%             else
                rk_to_rj = R1-R0;

                rjk = norm(rk_to_rj); %sqrt(sum((r(jj,:)-r(kk,:)).^2))
                rjk_hat = (rk_to_rj)/rjk;
                rjkrjk = rjk_hat*rjk_hat';
                %             keyboard
                if abs(rjk)>1e-15
                    Ajk = exp(1i*ko*rjk)/rjk*(ko^2*(rjkrjk - I) + (1i*ko*rjk-1)/rjk^2*(3*rjkrjk - I)); %Draine & Flatau
                    Toep(i,j,k,1) = Ajk(1,1);
                    Toep(i,j,k,2) = Ajk(1,2);
                    Toep(i,j,k,3) = Ajk(1,3);
                    Toep(i,j,k,4) = Ajk(2,2);
                    Toep(i,j,k,5) = Ajk(2,3);
                    Toep(i,j,k,6) = Ajk(3,3);
                end

%             end


        end
    end
end

tOpAss = toc(tIni);

opCirculant = circulant_nop_const(-Toep);
op_out = fft_operator(opCirculant);

s_opCirc = whos('opCirculant');

fACPU = @(J)mv_AN_const(J, op_out, Mr, Mc, 1/alpha_LDR, 'notransp', idxS3, 0);

%% Incident field with exp(ik\omega) dependecy as in Draine and Flatau
krx = kvec(1).*r(:,:,:,1);
kry = kvec(2).*r(:,:,:,2);
krz = kvec(3).*r(:,:,:,3);

kr = krx+kry+krz;

expKr=exp(1i*kr);

EINC(:,:,:,1) = Eo(1).*expKr;
EINC(:,:,:,2) = Eo(2).*expKr;
EINC(:,:,:,3) = Eo(3).*expKr;


%%
clear Vrhs
Vrhs = EINC(idxS3);


tini2=tic;
tic
% [circ_N,circ_L_opToep] = level_1_parallel_func_N_neaten_test(-Toep,L,...
%     M,N,1/alpha_LDR,0,'off');
[circ_N,circ_L_opToep] = level_1_parallel_func_N(-Toep,L,...
    M,N,1/alpha_LDR,0,'off');
disp('Circulant approximation of N operator');
toc

% [circ_M_opToep,circ_2_N] = circ_2_level_fast_neaten(circ_L_opToep,L,M,N);
[circ_M_opToep,circ_2_N] = circ_2_level_fast_test(circ_L_opToep,L,M,N);
%
parfor i=1:L
    for j=1:M
        circ_2_inv_T{i,j} = inv(1/alpha_LDR*eye(3*N)-circ_2_N{i,j});
    end
end
disp('2-level assembly')
t_2level = toc(tini2)

s_circ_2 = whos('circ_2_inv_T');

prec2_T = @(J) chan_2_mvp_for_parallel_idx(circ_2_inv_T,J,L,M,N,idxS3);

%% FFT accelerated
tol = 1e-5;
% % Solve without preconditioner
tiniD0=tic;
[vsolD0,~,~,~,resvecD0_GM] = pgmres(@(J)fACPU(J), Vrhs,2000, tol,1);%,@(J)prec2_T(J));
% [vsolD0,~,~,~,resvecD0] = bicgstab(@(J)fACPU(J), Vrhs, tol, 2000);
tendD0_GM=toc(tiniD0);
nIts_0 = length(resvecD0_GM);
fprintf('GMRES no preconditioner. Solve time = %.2f [sec] \n',tendD0_GM)
fprintf('Iteration count = %d \n',nIts_0);
%
tiniD0_BI=tic;
% [vsolD0,~,~,~,resvecD0] = pgmres(@(J)fACPU(J), Vrhs,2000, tol);%,@(J)prec2_T(J));
[vsolD0,~,~,~,resvecD0_BI] = bicgstab(@(J)fACPU(J), Vrhs, tol, 2000);
tendD0_BI=toc(tiniD0_BI);
nIts_0 = length(resvecD0_BI);
fprintf('BICG no preconditioner. Solve time = %.2f [sec] \n',tendD0_BI)
fprintf('Iteration count = %d \n',nIts_0);


% % % Solve with preconditioner
% tiniD1=tic;
% % [vsolD1,~,~,~,resvecD1] = pgmres(@(J)fACPU(J), Vrhs,2000, tol, 1,@(J)prec2_T(J));
% [vsolD1,~,~,~,resvecD1] = bicgstab(@(J)fACPU(J), Vrhs, tol, 2000,@(J)prec2_T(J));
% tendD1=toc(tiniD1);
% nIts_0 = length(resvecD1);
% fprintf('BICG with preconditioner. Solve time = %.2f [sec] \n',tendD1)
% fprintf('Iteration count = %d \n',nIts_0);
% % keyboard
% %
% % Solve with preconditioner
% tiniD1=tic;
% [vsolD1_GM,~,~,~,resvecD1_GM] = pgmres(@(J)fACPU(J), Vrhs,2000, tol, 1,@(J)prec2_T(J));
% % [vsolD1,~,~,~,resvecD1] = bicgstab(@(J)fACPU(J), Vrhs, tol, 2000,@(J)prec2_T(J));
% tendD1_GM=toc(tiniD1);
% nIts_0 = length(resvecD1_GM);
% fprintf('GMRES with preconditioner. Solve time = %.2f [sec] \n',tendD1_GM)
% fprintf('Iteration count = %d \n',nIts_0);

memOp = s_opCirc.bytes;
% memPrec = s_circ_2.bytes;

fprintf('Memory for operator = %.4e \n', memOp)
% fprintf('Memory for preconditioner = %.4e \n', memPrec)

fprintf('Num voxels = %.4e \n', numVox)

Chi = (refInd.^2-1)/(4*pi);

% vDDA = vsolD1./(Gram*Chi);
% vVIE = conj(vsolV1./(1j*omega*eo*(e_r-1)));

filename=['results_hex_DDA/hex_plate_fifth_DDA_no_prec_both_GMRES_BICG_refInd',REF_NAME{i_ref},'_sizeParam',num2str(sizeParam),...
    '_nPerLam',num2str(nPerLam)];
save(filename,'resvecD0_GM','resvecD0_BI','numVox','numVoxIce',...
    'tendD0_GM','tendD0_BI','tOpAss','s_opCirc')

% keyboard
% save(filename,'resvecD0_GM','numVox','numVoxIce',...
%     'tendD0_GM','tOpAss','s_opCirc')
% save(filename,'resvecD0','resvecD1','numVox','numVoxIce',...
%     'tendD0','tendD1','t_2level','tOpAss','s_opCirc','s_circ_2')
% save(filename,'resvecD1_GM','resvecD1_GM','numVox','numVoxIce',...
%     'tendD1_GM','tendD1_GM','t_2level','tOpAss','s_opCirc','s_circ_2')
% save(filename,'resvecD1_GM','resvecD1','numVox','numVoxIce',...
%     'tendD1_GM','tendD1','t_2level','tOpAss','s_opCirc','s_circ_2')
% save(filename,'resvecD1','numVox','numVoxIce',...
%     'tendD1','t_2level','tOpAss','s_opCirc','s_circ_2')

% filename=['results_hex/hex_plate_VIE_DDA_gmres_refInd','1p6','_sizeParam',num2str(sizeParam),...
%     '_nPerLam',num2str(nPerLam)];
% save(filename,'resvecV0','resvecV1','resvecD0','resvecD1','numVox','numVoxIce',...
%     'tendD0','tendD1','tendV0','tendV1')

clearvars -except SIZE_PARAM REF_RE i_param i_ref REF_NAME

    end
end
