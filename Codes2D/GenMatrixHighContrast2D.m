function GenMatrixHighContrast2D(filename, h0, p0, nmax, nz, offset)

% Purpose: custom script to generate stiffness matrices and partitioning
%          for the High contrast Poisson problem

% control the random seed for segmentation
rng(0);

Globals2D;

N=p0;

% choose Mesh/generation method
[VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h0, h0);
Nv = max(EToV(:));

StartUp2D;

BuildBCMaps2D;

% Build the create nz distinct zones with different diffusion parameters
mu = ones(K,1); mu1 = 1; mu2 = 1E-1;
%x_vals = linspace(-1,1,nz+1);
% introduce shift by one cell
x_vals = linspace(-1,1,nz+1) + offset/h0; x_vals(1) = -1; x_vals(end) = 1;
mu_vals = mod(1:nz,2)*mu1 + mod(2:nz+1,2)*mu2;
elements = (1:K)';
for i=1:nz
  sel = find( x_vals(i) <= ( VX(EToV(elements, 1)) + VX(EToV(elements, 2)) + VX(EToV(elements, 3)) )/3 & ( VX(EToV(elements, 1)) + VX(EToV(elements, 2)) + VX(EToV(elements, 3)) )/3 < x_vals(i+1))';
  mu(sel) = mu_vals(i)*mu(sel);
end

% Initialize solver and construct grid and metric
[A,M] = HighContrastIPDG2D(mu);

% set up Dirichlet boundary conditions
uD = zeros(Nfp*Nfaces, K);

% set up Neumann boundary conditions
% it's wrong and doesn't really matter
qN = zeros(Nfp*Nfaces, K);
%muN = ones(Nfp*Nfaces, K);


% evaluate boundary condition contribution to rhs
Aqbc = HighContrastIPDGbc2D(uD, qN, mu);

% set up right hand side forcing and adjust according to the jump in a
f = ones(Np,K);
rhs = MassMatrix*(J.*f) + Aqbc;

% % solve system
u = A\rhs(:);
u = reshape(u, Np, K);

% generate the elimination tree based on the coefficients
%elim_tree = GenElimTreeCoeff2D(4, mu, 'a');

% For plotting use
mu = ones(Np,1)*mu';
SurfDG2D(mu)

% generate eleimination tree
elim_tree = GenElimTree2D(nmax, 'PartitionMode', 'alternating');
%   elim_tree = GenElimTree2D(nmax, 'vertical');
% elseif strcmp(pmode, 'perpendicular')
%   elim_tree = GenElimTree2D(nmax, 'horizontal');
% else
%   elim_tree = GenElimTree2D(nmax, 'alternating');
% end
% PlotSparsityReorder(A, sep_tree);
%perm = ReorderElimTree(elim_tree, 'postorder');
% A = A(perm, perm);
% [L, U] = lu(A);
% figure
% spy(L);
% figure
% spy(A);
% convert elimination tree into the sequential format
elim_tree = ElimTree2SequentialTree(elim_tree);

% visualize the partitioning
%sel = elim_tree.inter(:,elim_tree.lsons == -1); sel = setdiff(unique(sel(:)), [0]);
%sel = elim_tree.bound(:,1); sel = setdiff(unique(sel(:)), [0]);
%v = zeros(size(rhs(:)));
%v(sel) = 1;
%v = reshape(v, Np, K);
%SurfDG2D(v)



% Save matrix, RHS and sep_tree in binary file
% Matlab stores matrices in CSC format
% If we want CSR format we can transpose before saving
% A = A'; transposed = true;
transposed = false;
b = rhs(:);
save(filename,'A','b','elim_tree','transposed','-v7.3','-nocompression');

return
