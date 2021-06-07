function GenMatrixHelmholtzCG2D(filename, h0, p0, k, nmax)

% custom script for Poisson Parameter study
Globals2D;
GlobalsCG2D;

N=p0;

% choose Mesh/generation method
[VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h0, h0);
%[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
%[Nv, VX, VY, K, EToV, BCType] = GenStructuredSquareQuadMesh2D(h0, h0);
Nv = max(EToV(:));

StartUp2D;
BuildBCMaps2D;
BuildGlobalMaps2D;

% Initialize solver and construct grid and metric
[A,M] = PoissonCG2D();
H = - k^2 * M + A;

% compute exact solution for Dirichlet boundary conditions
[uExact, fExact] = problem1(xCG,yCG);

% set up rhs
rhs = M(1:nFree, 1:nTotal)*fExact - H(1:nFree, nFree+1:nTotal)*uExact(nFree+1:nTotal);
S = A(1:nFree, 1:nFree);
M = M(1:nFree, 1:nFree);
u = H(1:nFree, 1:nFree) \ rhs(:);
uCG = [u; uExact(nFree+1:nTotal)];
SurfCG2D(uCG)
%uCG = [u; uExact(nFree+1:nTotal)];


elim_tree = GenElimTreeCG2D(nmax, 'PartitionMode', 'alternating');
%[perm, iperm] = ElimTree2Perm(elim_tree, 'postorder');
%figure
%spy(A(perm,perm));
%VisualizeTreeCG2D(elim_tree);
elim_tree = ElimTree2SequentialTree(elim_tree);


% Save matrix, RHS and sep_tree in binary file
% Matlab stores matrices in CSC format
% If we want CSR format we can transpose before saving
% A = A'; transposed = true;
A = H(1:nFree, 1:nFree);
transposed = false;
b = rhs(:);
save(filename,'A','S','M','b','elim_tree','transposed','-v7.3','-nocompression');

return
end

% define analytical solutions
function [u, f] = problem1(x,y)
  u = zeros(size(x));
  f = ones(size(x));
end