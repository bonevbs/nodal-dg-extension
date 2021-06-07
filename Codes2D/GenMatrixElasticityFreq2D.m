function GenMatrixElasticity2D(filename, h0, p0, k, icase, nmax)

% Purpose: custom script to generate stiffness matrices and partitioning
%          for the Helmholtz problem

Globals2D;

N=p0;

% parse input
if strcmp(icase, 'case 1') || strcmp(icase, 'heterogeneous')
  [VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h0, h0, 'Dirichlet');
elseif strcmp(icase, 'case 2')
  [VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h0, h0, 'channel');
else
  error("Unknown geometry.")
  %[VX, VY, K, EToV, BCType] = GenCircleMesh2D(h0);
end

%Nv = max(EToV(:));

StartUp2D;

BuildBCMaps2D;

% set up Dirichlet boundary conditions
uxD = zeros(Nfp*Nfaces, K);
uyD = zeros(Nfp*Nfaces, K);
% set up Neumann boundary conditions
qxN = zeros(Nfp*Nfaces, K);
qyN = zeros(Nfp*Nfaces, K);

% set problem specific settings
switch icase
  case 'case 1'
    uxD(mapD) = 0;
    uyD(mapD) = 0;
    qxN(mapN) = 0;
    qyN(mapN) = 1;
    % set parameters
    mu = ones(Np,K);
    lambda = ones(Np,K);
    fx = ones(Np,K);
    fy = zeros(Np,K);
    
  case 'case 2'
    uxD(mapD) = 0;
    uyD(mapD) = 0;
    qxN(mapN) = 0;
    qyN(mapN) = 1;
    % set parameters
    mu = ones(Np,K);
    lambda = ones(Np,K);
    fx = ones(Np,K);
    fy = zeros(Np,K);
  case 'heterogeneous'
    [mu, lambda] = heterogeneous(2, [1;0], 0, 'alternating 3', 10, 10);
    fx = zeros(size(x));
    fy = ones(size(x));
end

% Initialize solver and construct grid and metric
[Axx,Axy,Ayx,Ayy,M] = ElasticityIPDG2D(mu, lambda);
A = [Axx, Axy; Ayx, Ayy];
M = blkdiag(M,M);
H = - k^2 * M + A;

% evaluate boundary condition contribution to rhs
[bcx, bcy] = ElasticityIPDGbc2D(uxD,uyD,qxN,qyN,mu,lambda);

% set up right hand side forcing and adjust according to the jump in a
fx = - MassMatrix*(J.*fx); fy = - MassMatrix*(J.*fy);
rhs = [fx(:); fy(:)] + [bcx(:); bcy(:)];

% % solve system
u = H\rhs(:);
ux = u(1:Np*K); uy = u(Np*K+1:end);
ux = reshape(ux, Np, K); uy = reshape(uy, Np, K);

% For plotting use
%a = sqrt(ux.^2 + uy.^2);
%theta = atan2(uy,ux);
%figure
%SurfDG2D(a)
%figure
%SurfDG2D(theta)
figure
SurfDG2D(ux)
figure
SurfDG2D(uy)
figure
SurfDG2D(mu)
figure
SurfDG2D(lambda)
%figure
%QuiverDG2D(ux, uy);

% generate eleimination tree
elim_tree = GenElimTree2D(nmax, 'PartitionMode', "alternating");
%VisualizeTree2D(elim_tree)
elim_tree = StackElimTree(elim_tree, Np*K, 2);
perm = ElimTree2Perm(elim_tree, 'postorder');
%figure
%spy(A(perm,perm));
elim_tree = ElimTree2SequentialTree(elim_tree);


% Save matrix, RHS and sep_tree in binary file
% Matlab stores matrices in CSC format
% If we want CSR format we can transpose before saving
% A = A'; transposed = true;
transposed = false;
b = rhs(:);
A = H;
save(filename,'A','M','b','elim_tree','transposed','-v7.3','-nocompression');

return
end

function [mu, lambda] = heterogeneous(nzones, n, offset, mode, cmu, clambda)

Globals2D;

mu = zeros(Np,K);
lambda = zeros(Np,K);

n = n/norm(n);
VP = VX*n(1) + VY*n(2);
cmin = min(VP, [], 'all');
cmax = max(VP, [], 'all');

sep_left = -1;
for izone = 1:nzones
  sep_right = (cmax-cmin)/nzones*izone - 1 + offset;
  sel = ( ( VP(EToV(:, 1)) + VP(EToV(:, 2)) + VP(EToV(:, 3)) )/3 > sep_left ) & ...
    ( ( VP(EToV(:, 1)) + VP(EToV(:, 2)) + VP(EToV(:, 3)) )/3 <= sep_right );
  if strcmp(mode, 'increasing')
    mu(:,sel) = cmu*izone;
    lambda(:,sel) = clambda*izone;
  elseif strcmp(mode, 'alternating 1')
    mu(:,sel) = mod(izone,2)*cmu + mod(izone+1,2)*1;
    lambda(:,sel) = mod(izone,2)*clambda + mod(izone+1,2)*1;
  elseif strcmp(mode, 'alternating 2')
    mu(:,sel) = mod(izone+1,2)*cmu + mod(izone,2)*1;
    lambda(:,sel) = mod(izone,2)*clambda + mod(izone+1,2)*1;
  elseif strcmp(mode, 'alternating 3')
    mu(:,sel) = rand*cmu;
    lambda(:,sel) = rand*clambda;
  else
    error('Unknown mode')
  end
  sep_left = sep_right;
end

end

function [mu, lambda] = heterogeneous_bad(nzones, n, offset, mode, cmu, clambda)

Globals2D;

mu = ones(Np,K);
lambda = ones(Np,K);

n = n/norm(n);
coord = x*n(1) + y*n(2);
cmin = min(coord, [], 'all');
cmax = max(coord, [], 'all');

sep_left = -Inf;
for izone = 0:nzones
  sep_right = (cmax-cmin)/nzones*izone - 1 + offset;
  sel = ( coord > sep_left ) & ( coord <= sep_right );
  if strcmp(mode, 'increasing')
    mu(sel) = cmu*izone;
    lambda(sel) = clambda*izone;
  elseif strcmp(mode, 'alternating 1')
    mu(sel) = mod(izone,2)*cmu + mod(izone+1,2)*1;
    lambda(sel) = mod(izone,2)*clambda + mod(izone+1,2)*1;
  elseif strcmp(mode, 'alternating 2')
    mu(sel) = mod(izone+1,2)*cmu + mod(izone,2)*1;
    lambda(sel) = mod(izone,2)*clambda + mod(izone+1,2)*1;
  elseif strcmp(mode, 'alternating 3')
    mu(sel) = rand*cmu;
    lambda(sel) = rand*clambda;
    if (izone == nzones)
      mu(sel) = 0;
    end
  else
    error('Unknown mode')
  end
  sep_left = sep_right;
end

end