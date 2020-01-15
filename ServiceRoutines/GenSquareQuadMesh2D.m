function [VX, VY, K, EToV, BC] = GenSquareQuadMesh2D(Nx, Ny)

% function [VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h0)
% Purpose  : Generate 2D mesh on a square domain;   
% written by Boris Bonev

X = (-1:1/Nx:1);
Y = (-1:1/Ny:1);

[VX, VY] = meshgrid(X,Y);
Nv = length(VX);

EToV = delaunay(VX,VY);
K = size(EToV,1);

% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);

% flatten the coordinate vectors
VX = VX(:)';  VY = VY(:)';

% Boundary conditions
[EToE,EToF]= tiConnect2D(EToV);

% Dirichlet BC flag
Dirichlet = 6;

% initialize BC flags
BC = zeros(K,3);

% loop over elements
for i=1:K
    
    % loop over faces
    for j=1:3
        
        % if no neighbor
        if (EToE(i,j) == i)
            BC(i,j) = Dirichlet;
        end
    end
end

return
