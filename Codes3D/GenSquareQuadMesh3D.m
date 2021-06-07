function [VX, VY, VZ, K, EToV] = GenSquareQuadMesh3D(Nx, Ny, Nz)

% function [VX, VY, VZ, K, EToV] = GenSquareQuadMesh3D(Nx, Ny, Nz)
% Purpose  : Generate 3D mesh on a square domain;   
% written by Boris Bonev

X = (0:1/Nx:1);
Y = (0:1/Ny:1);
Z = (0:1/Nz:1);

[VX, VY, VZ] = meshgrid(X,Y,Z);
Nv = length(VX);

% flatten the coordinate vectors
VX = VX(:)';  VY = VY(:)'; VZ = VZ(:)';

EToV = delaunay(VX,VY,VZ);
K = size(EToV,1);

% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));  az = VZ(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));	 bz = VZ(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));	 cz = VZ(EToV(:,3));
dx = VX(EToV(:,4)); dy = VY(EToV(:,4));	 dz = VZ(EToV(:,4));

Dx = (cy-ay).*(dz-az) - (cz-az).*(dy-ay);
Dy = (cz-az).*(dx-ax) - (cx-ax).*(dz-az);
Dz = (cx-ax).*(dy-ay) - (cy-ay).*(dx-ax);

D = (bx-ax).*Dx + (by-ay).*Dy + (bz-az).*Dz;

i = find(D<0);
EToV(i,:) = EToV(i,[1 2 4 3]);


return
