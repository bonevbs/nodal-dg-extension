function GenMatrixPoisson3D(filename, h0, p0, nmax)
% custom script for Poisson Parameter study
Globals3D;

N=p0;

% choose Mesh/generation method
%[VX, VY, K, EToV] = GenSquareQuadMesh2D(h0, h0);
%[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
[VX, VY, VZ, K, EToV] = GenSquareQuadMesh3D(h0, h0, h0);
Nv = max(EToV(:));

StartUp3D;

%PlotMesh2D();

% Build boundary conditions when using unstructured mehs
BCType = (EToE == (1:K)' * [1 1 1 1]);
BCType = BCType.*Dirichlet;

BuildBCMaps3D;

% Initialize solver and construct grid and metric
[A,M] = PoissonIPDG3D();

% Plot mesh
% PlotMesh3D();

% set up boundary condition 
xbc = Fx(mapB); ybc = Fy(mapB); zbc = Fz(mapB);
ubc = zeros(Nfp*Nfaces*K,1);
ubc(mapB) = sin(pi*xbc).*sin(pi*ybc).*sin(pi*zbc);

% form right hand side contribution from boundary condition
Abc = PoissonIPDGbc3D(ubc);

% evaluate forcing function
f = -3*(pi^2)*sin(pi*x).*sin(pi*y).*sin(pi*z); 

% set up right hand side for variational Poisson equation
rhs = M*(-f(:)) + Abc(:);

% % solve system
%u = A\rhs(:);
%u = reshape(u, Np, K);

% % For plotting use
%[TRI,xout,yout,uout,interp] = PlotField2D(2*N,x,y,u);
%clf;
%PlotDomain2D();
%hold on;
%PlotContour2D(TRI,xout,yout,uout,linspace(0,-12000,12));
%hold off;

elim_tree = GenElimTree3D(nmax, 'PartitionMode', 'alternating');
[perm, iperm] = ElimTree2Perm(elim_tree, 'postorder');
elim_tree = ElimTree2SequentialTree(elim_tree);

% spy(A)
% figure
% spy(A(perm, perm))


% Save matrix, RHS and sep_tree in binary file
% We make use of the fact that Matlab stores matrices in CSC format
% So matrix is stored in CSR format as we want it
transposed = false;
b = rhs(:); x = x(:); y = y(:); z = z(:);
save(filename,'A','b','elim_tree','x','y','z','transposed','-v7.3','-nocompression');

return