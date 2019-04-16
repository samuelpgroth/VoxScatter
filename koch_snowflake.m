% clear

%% Create Koch snowflake geometry
% I pinched this first little piece of code (generating the vertices) from
% https://people.cs.clemson.edu/~goddard/handouts/math360/notes11.pdf

P = [ 0 0; 1 0; cos(-pi/3), sin(-pi/3); 0 0  ];   % a triangle
for iteration=1:5
    newP = zeros( size(P,1)*4+1, 2);
    for i=1:size(P,1)-1
        newP(4*i+1,:) = P(i,:);
        newP(4*i+2,:) = (2*P(i,:) + P(i+1,:) )/3;
        link = P(i+1,:)-P(i,:);
        ang = atan2( link(2), link(1) );
        linkLeng = sqrt( sum(link.^2) );
        newP(4*i+3,:) = newP(4*i+2,:) + (linkLeng/3)*[ cos(ang+pi/3), sin(ang+pi/3) ];
        newP(4*i+4,:) = (P(i,:) + 2*P(i+1,:) )/3;
    end
    newP( 4*size(P,1)+1,:) = P(size(P,1),:);
    P = newP;
end

%% Now I go about creating a bounding box and the voxelized grid. 

% Establish bounding box around snowflake
xmin = min(P(:,1));
ymin = min(P(:,2));
xmax = max(P(:,1));
ymax = max(P(:,2));

dom_y = ymax - ymin;

% Shift so centre is at the origin and then scale so that radius is 1
P(:,1) = P(:,1) - (xmax+xmin)./2;
P(:,2) = P(:,2) - (ymax+ymin)./2;

P = P./ (dom_y/2);

xmin = min(P(:,1));
ymin = min(P(:,2));
xmax = max(P(:,1));
ymax = max(P(:,2));

dom_x = xmax - xmin;
dom_y = ymax - ymin;
dom_z = dom_y * aspectRatio;

%% Create uniform grid of bounding box, including endpoints

% Preferential dimension: 
%this is the direction in which we ensure that the voxels match the
%geometry
h_pref = dom_y; % here it's the width (y-direction)
% figure out resolution
res_temp = lambda_int/nPerLam; % provisional resolution
N = ceil(h_pref/res_temp);
res = h_pref/N;

% Generate voxels
[r] = generatedomain_new(res,dom_x/2,dom_y/2,dom_z/2); 
[L,M,N,~]=size(r);

% Nx = 100; % number of voxels in x-direction
% h  = dx/Nx; % voxel dimension
% Ny = ceil(dy_temp/h);
% dy = Ny*h;
% Nz = 10;
% 
% % Readjust ymin and ymax, given new dy
% ymin = ymin - (dy-dy_temp)/2;
% ymax = ymax + (dy-dy_temp)/2;
% zmin = -Nz/2 * h;
% zmax = Nz/2 * h;
% 
% 
% x = linspace(xmin+h/2,xmax-h/2,Nx);
% y = linspace(ymin+h/2,ymax-h/2,Ny);
% z = linspace(zmin+h/2,zmax-h/2,Nz);
% 
% for i = 1:Nx
%     for j = 1:Ny
%         for k = 1:Nz
%             r(i,j,k,1) = x(i);
%             r(i,j,k,2) = y(j);
%             r(i,j,k,3) = z(k);
%         end
%     end
% end

% plot geometry
xd = r(:,:,:,1);
yd = r(:,:,:,2);
zd = r(:,:,:,3);

% figure(1);
% plot3(xd(:), yd(:), zd(:), 's');
% axis equal;
% grid on;

% Inpolygon decides which of my voxels lie inside and outside the snowflake
[in,on] = inpolygon(xd(:),yd(:),P(:,1),P(:,2));

% figure
% plot3(xd(in),yd(in),zd(in),'r+') % points inside
% hold on
% plot3(xd(~in),yd(~in),zd(~in),'bo') % points outside
% axis equal

figure
plot3(xd(in),yd(in),zd(in),'.')
axis image


% %% For Karina's code where she just wants the dipoles inside the shape  
% xx = xd(in);
% yy = yd(in);
% zz = zd(in);
% r = [xx,yy,zz];



