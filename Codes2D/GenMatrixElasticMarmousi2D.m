function GenMatrixElasticMarmousi2D(filename, h0, p0, omega, mesh_type, nmax)

% Purpose: custom script to generate stiffness matrices and partitioning
%          for the Helmholtz problem

Globals2D;

N=p0;

if strcmp(mesh_type, 'regular')
  % load Marmousi data
  load('/Users/boris/Code/Matlab/nodal-dg-extension/Grid2D/marmousi.mat');
  % determine dimensions
  xmin = min(xM, [], 'all');
  xmax = max(xM, [], 'all');
  ymin = min(yM, [], 'all');
  ymax = max(yM, [], 'all');
  ymax = -450; % cut off the water
  h0 = (ymax - ymin)/h0;
  Nx = round_even((xmax - xmin)/h0);
  Ny = round_even((ymax - ymin)/h0);
  % use rectangular mesh
  [VX, VY, K, EToV, BCType] = GenRectangleMesh2D(Nx, Ny, [xmin,ymin;xmax,ymax], 'Dirichlet');
elseif strcmp(mesh_type, 'conforming')
  [VX, VY, K, EToV, BCType, Rho, Mu, Lambda] = GenMarmousiMesh2D(h0, 'bottom');
else
  error('Unknown mesh type.')
end

StartUp2D;

BuildBCMaps2D;

if strcmp(mesh_type, 'regular')
  Rho = 100^3 * 0.001 * interp2(XM,YM,rho,x,y);
  Vp = interp2(XM,YM,v_p,x,y);
  Vs = interp2(XM,YM,v_s,x,y);
  Mu = Rho.*Vs.^2;
  Lambda = Rho.*Vp.^2 - 2*Mu;
else
  Rho = ones(Np,1)*(Rho(:))';
Mu = ones(Np,1)*(Mu(:))';
Lambda = ones(Np,1)*(Lambda(:))';
end
figure
SurfDG2D(Rho)
figure
SurfDG2D(Mu)
figure
SurfDG2D(Lambda)

% set up Dirichlet boundary conditions
uxD = zeros(Nfp*Nfaces, K);
uyD = zeros(Nfp*Nfaces, K);
uxD(mapD) = 0;
uyD(mapD) = 0;
% set up Neumann boundary conditions
qxN = zeros(Nfp*Nfaces, K);
qyN = zeros(Nfp*Nfaces, K);

% set forcing term
Fx = zeros(size(x));
Fy = zeros(size(x));
xs = 8500;
ys = -1250;
rs = 100;
Fx = (x-xs).*exp(-((x - xs).^2 + (y - ys).^2)/(2*rs^2));
Fy = (y-ys).*exp(-((x - xs).^2 + (y - ys).^2)/(2*rs^2));
%Fx = sign(y - ys).*exp(-((x - xs).^2 + (y - ys).^2)/(2*rs^2));
%figure
%SurfDG2D(Fx)

% Initialize solver and construct grid and metric
[Axx,Axy,Ayx,Ayy,M] = ElasticityIPDG2D(Mu, Lambda);
A = [Axx, Axy; Ayx, Ayy];
Rho_d = sparse(1:(K*Np),1:(K*Np),Rho(:),K*Np,K*Np);
M = blkdiag(M*Rho_d,M*Rho_d);
%M = blkdiag(M,M);
H = - omega^2 * M + A;

% evaluate boundary condition contribution to rhs
[bcx, bcy] = ElasticityIPDGbc2D(uxD,uyD,qxN,qyN,Mu,Lambda);

% set up right hand side forcing and adjust according to the jump in a
Fx = - MassMatrix*(J.*Fx); Fy = - MassMatrix*(J.*Fy);
rhs = [Fx(:); Fy(:)] + [bcx(:); bcy(:)];

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
%figure
%QuiverDG2D(ux, uy);
% figure
% QuiverDG2D(ux, uy);
% figure
% ContfdDG2D(a, ux, uy, 1e-2);

% generate eleimination tree
elim_tree = GenElimTree2D(nmax, 'PartitionMode', "aspect");
figure
VisualizeTree2D(elim_tree)
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

function S = round_even(S)
% round to nearest even integer.
idx = mod(S,2)>1;
S = floor(S);
S(idx) = S(idx)+1;
end
