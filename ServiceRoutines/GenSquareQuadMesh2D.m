function [VX, VY, K, EToV, BC] = GenSquareQuadMesh2D(Nx, Ny, varargin)

% function [VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h0)
% Purpose  : Generate 2D mesh on a square domain;   
% written by Boris Bonev

BCType = 'Dirichlet'; % Use Dirichlet by default
if nargin > 2
    if strcmp(varargin{1}, 'Dirichlet') || strcmp(varargin{1}, 'Neumann')...
        || strcmp(varargin{1}, 'corner') || strcmp(varargin{1}, 'channel')
      BCType = varargin{1};
    else
      error('GenSquareQuadMesh2D: Unknown boundary type specified.')
    end
end

X = (-1:2/Nx:1);
Y = (-1:2/Ny:1);

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

% Boundary conditions
[EToE,EToF]= tiConnect2D(EToV);

% initialize BC flags
BC = zeros(K,3);

% loop over elements
for i=1:K
    
  % loop over faces
  for j=1:3
        
    % if no neighbor
    if (EToE(i,j) == i)
      if strcmp(BCType, 'Dirichlet')
        BC(i,j) = 6;
      elseif strcmp(BCType, 'Neumann')
        BC(i,j) = 7;
      elseif strcmp(BCType, 'corner')
        if any(VX(EToV(i,1:3)) == -1) || any(VY(EToV(i,1:3)) == -1)
          BC(i,j) = 6;
        else
          BC(i,j) = 7;
        end
      elseif strcmp(BCType, 'channel')
        if ( VX(EToV(i,j)) == -1 && VX(EToV(i,mod(j,3)+1)) == -1 ) ||...
          ( VX(EToV(i,j)) == 1 && VX(EToV(i,mod(j,3)+1)) == 1 )
          BC(i,j) = 6;
        elseif ( VY(EToV(i,j)) == -1 && VY(EToV(i,mod(j,3)+1)) == -1 ) ||...
          ( VY(EToV(i,j)) == 1 && VY(EToV(i,mod(j,3)+1)) == 1 )
          BC(i,j) = 7;
        else
          error("Node should not be on the boundary.")
        end
      end
    end
  end
end

% flatten the coordinate vectors
VX = VX(:)';  VY = VY(:)';

return
