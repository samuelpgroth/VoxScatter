%------------------------------------------------------%
%      TEST getOPERATORS() function                    %
%------------------------------------------------------%

% load RHBM.r  
load('body.mat');

% set desired frequency 
f = 3e8;

%choose appropriate solver for singular integrals 
singular = 'Directfn';
%singular = 'Demcem';

% select solver's structure
type = 'F';                             % T - Toeplitz; C - Circulant; 
                                        % FFT_circulant is default (any other symbol)
                                        
                                        

% compute N operator
[N] = getOPERATORS(r, f, 'N', type, singular);        % test N operator

% compute K operator 
[K] = getOPERATORS(r, f, 'K', type, singular);        % test K operator