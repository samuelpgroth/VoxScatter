function [Esca] = E_field_Nop_const_Sam(J,fN,Gram,freq,Einc,GPU_flag)
% Returns only scattered field, rather than total field
%%   Generation of E field 
% _________________________________________________________________________
%
%       Applies the N operator to a current density vector J
%       and adds the incident Electric field to generate
%       the total fields
%
% _________________________________________________________________________
%
%% INPUT
%   J:          Volumetric current density (LxMxNx3)
%   fN:         FFT-Circulant of N operator
%   Gram:       Gram matrix
%   freq:       frequency (in Hz) 
%   Einc:       incident field (LxMxNx3) -- point matching!
%   GPU_flag:   applies GPU if 1, no GPU if 0
%
%
%% OUTPUT
%   Etot: total field (LxMxNx3) -- point matching!
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________


% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 4 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if((nargin < 5) || isempty(Einc))
    Einc = 0*J; % zero incident field
end
if((nargin < 6) || isempty(GPU_flag))
    GPU_flag = 0; % no GPU
end

% Check if currents are zero for fast return
if (nnz(J) == 0)
    Etot = Einc;
    return
end


% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------

% EM related constants
mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
omega = 2*pi*freq;

% domain dimensions
[L, M, N, ~] = size(J);

% positions is all the domain
idx = 1:3*L*M*N;

% ---------------------------------------------------------------------
% Call the L operator
% ---------------------------------------------------------------------

% we pass J instead of r, since r in mv_L is only use to get dimensions
[Esca] = mv_L_const(J, fN, J, Gram, idx, GPU_flag);
Esca = reshape(Esca,L,M,N,3);

% -------------------------------------------------------------------------
% add incident field contribution to the radiated field (in Etot)
% -------------------------------------------------------------------------

% divide by j*omega*eo and scale by voxel volume to obtain point fields
Esca = Esca./(Gram*1j*omega*eo); % Etotal = Einc + (E - dV*J)/(dV*j*omega*eo)

