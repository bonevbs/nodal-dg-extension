function [VX, VY, K, EToV, BC] = GenRectangleMesh2D(Nx, Ny, dims, varargin)

% function [VX, VY, K, EToV, BCType] = GenRectangleMesh2D(h0)
% Purpose  : Generate 2D mesh on a rectangular domain;
%    dims  : Bounding box [xmin,ymin; xmax,ymax]
% written by Boris Bonev

BCType = 'Dirichlet'; % Use Dirichlet by default
if nargin > 3
    if strcmp(varargin{1}, 'Dirichlet') || strcmp(varargin{1}, 'Neumann')...
        || strcmp(varargin{1}, 'corner') || strcmp(varargin{1}, 'channel')...
        || strcmp(varargin{1}, 'bottom') || strcmp(varargin{1}, 'Dirac')
      BCType = varargin{1};
    else
      error('GenSquareQuadMesh2D: Unknown boundary type specified.')
    end
end

lx = dims(2,1) - dims(1,1);
ly = dims(2,2) - dims(1,2);

xmin = dims(1,1);
xmax = dims(2,1);
ymin = dims(1,2);
ymax = dims(2,2);

X = (xmin:lx/Nx:xmax);
Y = (ymin:ly/Ny:ymax);

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
        if any(VX(EToV(i,1:3)) <= xmin) || any(VY(EToV(i,1:3)) <= ymin)
          BC(i,j) = 6;
        else
          BC(i,j) = 7;
        end
      elseif strcmp(BCType, 'channel')
        if ( VX(EToV(i,j)) <= xmin && VX(EToV(i,mod(j,3)+1)) <= xmin ) ||...
          ( VX(EToV(i,j)) >= xmax && VX(EToV(i,mod(j,3)+1)) >= xmax )
          BC(i,j) = 6;
        elseif ( VY(EToV(i,j)) <= ymin && VY(EToV(i,mod(j,3)+1)) <= ymin ) ||...
          ( VY(EToV(i,j)) >= ymax && VY(EToV(i,mod(j,3)+1)) >= ymax )
          BC(i,j) = 7;
        else
          error("Node should not be on the boundary.")
        end
      elseif strcmp(BCType, 'bottom')
        if ( VX(EToV(i,j)) <= xmin && VX(EToV(i,mod(j,3)+1)) <= xmin ) ||...
          ( VX(EToV(i,j)) >= xmax && VX(EToV(i,mod(j,3)+1)) >= xmax )
          BC(i,j) = 7;
        elseif ( VY(EToV(i,j)) >= ymax && VY(EToV(i,mod(j,3)+1)) >= ymax )
          BC(i,j) = 7;
        elseif ( VY(EToV(i,j)) <= ymin && VY(EToV(i,mod(j,3)+1)) <= ymin )
          BC(i,j) = 6;
        else
          error("Node should not be on the boundary.")
        end
      elseif strcmp(BCType, 'Dirac')
        if ( VX(EToV(i,j)) == (xmax-xmin)/2) || ( VY(EToV(i,j)) == ymax && VY(EToV(i,mod(j,3)+1)) == ymax )
          BC(i,j) = 6;
        else
          BC(i,j) = 7;
        end
      end
    end
  end
end

% flatten the coordinate vectors
VX = VX(:)';  VY = VY(:)';

return
