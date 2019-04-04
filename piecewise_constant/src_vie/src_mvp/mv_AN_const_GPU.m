function [JOut] = mv_AN_const_GPU(JIn0, fN, Mr, Mc, Gram, idx)
%%   M-V product of operator A_N via FFT
% _________________________________________________________________________
%
%       Applies the AN operator to a current vector
%
% _________________________________________________________________________
%
%% INPUT
%   JIn:            Current 
%   fN:             FFT-Circulant of N operator
%   Mr:             relative permittivity or permeability 
%   Mc:             suceptibility
%   Gram:           Gram matrix 
%   transp_flag:    flag for choosing the m-v product
%                   'transp':   AN' * JIn0
%                   'notransp': AN * JIn0
%   idx:            index with local coordinates (non-air voxels)
%   GPU_flag:       applies GPU if 1, no GPU if 0
%
%
%% OUTPUT
%   JOut:           Current
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
% Prepare data
% -------------------------------------------------------------------------

% fft dimensions
[LfN, MfN, NfN, ~] = size(fN);

% domain dimensions
[L, M, N] = size(Mr);


% allocate space
JOut = gpuArray.zeros(L, M, N, 3);
JIn = gpuArray.zeros(L, M, N, 3);
% send to gpu, and translate from local (idx) to global (L,M,N) coordinates
%JIn(idx) = gpuArray(JIn0);
JIn(idx) = JIn0;



    
    % ---------------------------------------------------------------------
    % apply fft and mv-op for each of the components of JIn
    % ---------------------------------------------------------------------
       
    % x component of JIn, store contribution on 3 components of Jout
    fJ = fftn(JIn(:,:,:,1),[LfN, MfN, NfN]);
    Jout1 = fN(:,:,:,1) .* fJ;
    Jout2 = fN(:,:,:,2) .* fJ;
    Jout3 = fN(:,:,:,3) .* fJ;

    % y component of JIn, add contribution on 3 components of Jout
    fJ = fftn(JIn(:,:,:,2),[LfN, MfN, NfN]);
    Jout1 = Jout1 + fN(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fN(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fN(:,:,:,5) .* fJ;
    
    % z component of JIn, add contribution on 3 components of Jout
    fJ = fftn(JIn(:,:,:,3),[LfN, MfN, NfN]);
    Jout1 = Jout1 + fN(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fN(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fN(:,:,:,6) .* fJ;
    
    % apply ifft, multiply by material properties and Gram
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = Gram .* Mr .* JIn(:,:,:,1) - Mc .* Jout1(1:L,1:M,1:N);
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = Gram .* Mr .* JIn(:,:,:,2) - Mc .* Jout2(1:L,1:M,1:N);
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = Gram .* Mr .* JIn(:,:,:,3) - Mc .* Jout3(1:L,1:M,1:N);


% -------------------------------------------------------------------------
% Return local coordinates related to material positions
% -------------------------------------------------------------------------

% get from GPU
%JOut = gather(JOut(idx));
JOut = JOut(idx);
% clear gpu data
clear JIn; clear Jout1; clear Jout2; clear Jout3; clear fJ;


