function GenMatrixHelmholtz2D(filename, h0, p0, k, geometry, nmax)

% Purpose: custom script to generate stiffness matrices and partitioning
%          for the Helmholtz problem

Globals2D;

N=p0;

% choose Mesh/generation method
if strcmp(geometry, 'square')
  [VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h0, h0, 'Dirichlet');
elseif strcmp(geometry, 'circle')
  [VX, VY, K, EToV, BCType] = GenCircleMesh2D(h0, 'Dirichlet');
elseif strcmp(geometry, 'guitar')
  [VX, VY, K, EToV, BCType] = GenGuitarMesh2D(h0, 'mixed');
else
  error("Unknown geometry.")
end

Nv = max(EToV(:));

StartUp2D;

BuildBCMaps2D;

% Initialize solver and construct grid and metric
[A,M] = PoissonIPDG2D();
H = - k^2 * M + A;

% set up Dirichlet boundary conditions
uD = zeros(Nfp*Nfaces, K);

% set up Neumann boundary conditions
% it's wrong and doesn't really matter
qN = zeros(Nfp*Nfaces, K);
%muN = ones(Nfp*Nfaces, K);


% evaluate boundary condition contribution to rhs
Aqbc = PoissonIPDGbc2D(uD, qN);

% set up right hand side forcing and adjust according to the jump in a
f = ones(Np,K); % constant one
%f = (x == 0 & y == 0); % Dirac delta
%f = 1*exp(- (x.^2 + y.^2)/(2*0.0000001));
rhs = MassMatrix*(J.*f) + Aqbc;

% % solve system
u = H\rhs(:);
u = reshape(u, Np, K);

figure
% For plotting use
SurfDG2D(u)

% generate eleimination tree
elim_tree = GenElimTree2D(nmax, 'PartitionMode', 'aspect');
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
%VisualizeTree2D(elim_tree)
elim_tree = ElimTree2SequentialTree(elim_tree);


% Save matrix, RHS and sep_tree in binary file
% Matlab stores matrices in CSC format
% If we want CSR format we can transpose before saving
% A = A'; transposed = true;
transposed = false;
b = rhs(:);
A = H;
save(filename,'A','M','b','elim_tree','k','transposed','-v7.3','-nocompression');

return
