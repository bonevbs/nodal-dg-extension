function GenMatrixElasticity2D(filename, h0, p0, icase, nmax)

% Purpose: custom script to generate stiffness matrices and partitioning
%          for the Helmholtz problem

Globals2D;

N=p0;

% parse input
if strcmp(icase, 'forced 1') || strcmp(icase, 'heterogeneous') || strcmp(icase, 'bimaterial bar')
	[VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h0, h0, 'Dirichlet');
elseif strcmp(icase, 'forced 2')
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

switch icase
  case 'forced 1'
    [uxA, uyA, fx, fy, mu, lambda] = solution1(x,y);
    uxD(mapD) = uxA(vmapD);
    uyD(mapD) = uyA(vmapD);
  case 'forced 2'  
    [uxA, uyA, qxA, qyA, fx, fy, mu, lambda] = solution2(x,y);
    uxD(mapD) = uxA(vmapD);
    uyD(mapD) = uyA(vmapD);
    qxN(mapN) = qxA(vmapN);
    qyN(mapN) = qyA(vmapN);
  case 'heterogeneous'
    [mu, lambda] = heterogeneous(2, 'vertical');
    fx = zeros(size(x));
    fy = ones(size(x));
  case 'bimaterial bar'
    e1 = 1; e2 = 10;
    mu1 = e1/2; mu2 = e2/2;
    % interface is at 0
    alpha = e2/(e2 + e1);
    sel2 = ( VX(EToV(:, 1)) + VX(EToV(:, 2)) + VX(EToV(:, 3)) )/3 > 0;
    mu = mu1*ones(size(x));
    lambda = zeros(size(x));
    mu(:,sel2) = mu2;
    uxA = (1+x)*alpha;
    uxA(:,sel2) = 1+e1/e2*(x(:,sel2) - 1)*alpha;
    uyA = zeros(size(x));
    %lambda(:,sel2) = 10;
    fx = zeros(size(x));
    fy = zeros(size(x));
    uxD(mapD) = uxA(vmapD);
    uyD(mapD) = uyA(vmapD);
end

% Initialize solver and construct grid and metric
[Axx,Axy,Ayx,Ayy,M] = ElasticityIPDG2D(mu, lambda);
A = [Axx, Axy; Ayx, Ayy];
M = blkdiag(M,M);

% evaluate boundary condition contribution to rhs
[bcx, bcy] = ElasticityIPDGbc2D(uxD,uyD,qxN,qyN,mu,lambda);

% set up right hand side forcing and adjust according to the jump in a
fx = - MassMatrix*(J.*fx); fy = - MassMatrix*(J.*fy);
rhs = [fx(:); fy(:)] + [bcx(:); bcy(:)];

% % solve system
u = A\rhs(:);
ux = u(1:Np*K); uy = u(Np*K+1:end);
ux = reshape(ux, Np, K); uy = reshape(uy, Np, K);

% For plotting use
a = sqrt(ux.^2 + uy.^2);
%theta = atan2(uy,ux);
%figure
%SurfDG2D(a)
%figure
%SurfDG2D(theta)
figure
SurfDG2D(ux)
figure
SurfDG2D(uy)
%figure
%SurfDG2D(mu)
%figure
%SurfDG2D(lambda)
%figure
%QuiverDG2D(ux, uy);
% figure
% QuiverDG2D(ux, uy);
% figure
% ContfdDG2D(a, ux, uy, 1e-2);

% generate eleimination tree
elim_tree = GenElimTree2D(nmax, 'PartitionMode', "alternating");
%VisualizeTree2D(elim_tree)
elim_tree = StackElimTree(elim_tree, Np*K, 2);
%perm = ElimTree2Perm(elim_tree, 'postorder');
%figure
%spy(A(perm,perm));
elim_tree = ElimTree2SequentialTree(elim_tree);


% Save matrix, RHS and sep_tree in binary file
% Matlab stores matrices in CSC format
% If we want CSR format we can transpose before saving
% A = A'; transposed = true;
transposed = false;
b = rhs(:);
save(filename,'A','M','b','elim_tree','transposed','-v7.3','-nocompression');

return
end

% define analytical solutions
function [ux, uy, fx, fy, mu, lambda] = solution1(x,y)
  c = 0;
  muc = 1; mux = 3; muy = 1;
  lambdac = 1; lambdax = 1; lambday = 2;
  ux = sin(pi*x).*sin(pi*y) + c;
  uy = sin(pi*x).*sin(pi*y) + c;
  mu = muc + mux*x + muy*y;
  lambda = lambdac + lambdax*x + lambday*y;
  fx = pi*lambdax*cos(pi*x).*sin(pi*y) + pi*lambdax*sin(pi*x).*cos(pi*y)...
    - pi^2*lambda.*sin(pi*x).*sin(pi*y) + pi^2*lambda.*cos(pi*x).*cos(pi*y)...
    + 2*pi*mux*cos(pi*x).*sin(pi*y) - 2*pi^2*mu.*sin(pi*x).*sin(pi*y)...
    + pi*muy*cos(pi*x).*sin(pi*y) + pi*muy*sin(pi*x).*cos(pi*y)...
    - pi^2*mu.*sin(pi*x).*sin(pi*y) + pi^2*mu.*cos(pi*x).*cos(pi*y);
  fy = pi*lambday*cos(pi*x).*sin(pi*y) + pi*lambday*sin(pi*x).*cos(pi*y)...
    - pi^2*lambda.*sin(pi*x).*sin(pi*y) + pi^2*lambda.*cos(pi*x).*cos(pi*y)...
    + 2*pi*muy*sin(pi*x).*cos(pi*y) - 2*pi^2*mu.*sin(pi*x).*sin(pi*y)...
    + pi*mux*cos(pi*x).*sin(pi*y) + pi*mux*sin(pi*x).*cos(pi*y)...
    - pi^2*mu.*sin(pi*x).*sin(pi*y) + pi^2*mu.*cos(pi*x).*cos(pi*y);
end

function [uxA, uyA, qxA, qyA, fx, fy, mu, lambda] = solution2(x,y)
    uxA = sin(pi*x).*cosh(pi*y);
    uyA = -cos(pi*x).*sinh(pi*y); 
    lambda = 1.0*ones(size(x));
    mu = 1.0*ones(size(x));
    fx = zeros(size(x));
    fy = zeros(size(x));
    qxA = 2*pi*sign(y).*mu.*sin(pi*x).*sinh(pi*y);
    qyA = -2*pi*sign(y).*mu.*cos(pi*x).*cosh(pi*y);
end

function [mu, lambda] = heterogeneous(nzones, orientation)

Globals2D;

mu = zeros(Np,K);
lambda = zeros(Np,K);

if strcmp(orientation, 'vertical')
  VP = VX;
elseif strcmp(orientation, 'horizontal')
  VP = VY;
else
  error('Unknown orientation')
end

sep_left = -1;
for izone = 1:nzones
  sep_right = 2/nzones*izone - 1;
  sel = ( ( VP(EToV(:, 1)) + VP(EToV(:, 2)) + VP(EToV(:, 3)) )/3 > sep_left ) & ...
    ( ( VP(EToV(:, 1)) + VP(EToV(:, 2)) + VP(EToV(:, 3)) )/3 <= sep_right );
  mu(:,sel) = izone;
  lambda(:,sel) = izone;
  sep_left = sep_right;
end

end
