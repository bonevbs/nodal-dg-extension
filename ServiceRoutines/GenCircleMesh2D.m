function [VX, VY, K, EToV, BC] = GenCircleMesh2D(h0, varargin)

% function [VX, VY, K, EToV] = MeshGenDistMesh2D()
% Purpose  : Generate 2D square mesh using DistMesh;   
% By Allan P. Engsig-Karup

% Parameters to set/define
%    fd     Distance function for mesh boundary
%    fh     Weighting function for distributing elements
%    h0     Characteristic length of elements
%    Bbox   Bounding box for mesh
%    param  Parameters to be used in function call with DistMesh

BCType = 'Dirichlet'; % Use Dirichlet by default
if nargin > 1
    if strcmp(varargin{1}, 'Dirichlet') || strcmp(varargin{1}, 'Neumann')...
        || strcmp(varargin{1}, 'mixed')
      BCType = varargin{1};
    else
      error('GenCircleMesh2D: Unknown boundary type specified.')
    end
end
if strcmp(BCType, 'Dirichlet')
  BCType = 6;
else
  BCType = 7;
end

fd = @(p) sqrt(sum(p.^2,2))-1;
fh = @huniform;
Bbox = [-1 -1; 1 1];
param = [];

% Call distmesh
[Vert,EToV]=distmesh2d(fd,fh,1/h0,Bbox,param);
VX = Vert(:,1); VY = Vert(:,2);
Nv = length(VX); K  = size(EToV,1);

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
            BC(i,j) = BCType;
        end
    end
end

% flatten the coordinate vectors
VX = VX(:)';  VY = VY(:)';
return
