function[r, idx, res, P, lambda_ext, lambda_int] = geometry(geom, refInd, sizeParam, nPerLam)

% Geometrical parameters (alter appropriate geometry)
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
    aspectRatio = 1/10; % ratio of column's height to it's width
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

%% Define some standard EM quantities
lambda_ext = 2*pi*a/sizeParam; % exterior wavelength of incident wave

%% DISCRETIZATION (into voxels)
% Preferential dimension: 
%this is the direction in which we ensure that the voxels match the
%geometry
h_pref = dom_x; % here it's the height (z-direction)
lambda_int = lambda_ext/real(refInd); % interior wavelength
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
