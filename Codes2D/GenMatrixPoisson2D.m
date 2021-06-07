function GenMatrixPoisson2D(filename, h0, p0, nmax)

% custom script for Poisson Parameter study
Globals2D;

N=p0;

% choose Mesh/generation method
[VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h0, h0);
%[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
%[Nv, VX, VY, K, EToV, BCType] = GenStructuredSquareQuadMesh2D(h0, h0);
Nv = max(EToV(:));

StartUp2D;

%PlotMesh2D();

BuildBCMaps2D;

% Initialize solver and construct grid and metric
[A,M] = PoissonIPDG2D();

% Plot mesh
% PlotMesh2D();

% set up Dirichlet boundary conditions
uD = zeros(Nfp*Nfaces, K);
[uD(mapD), ~] = solution1(Fx(mapD), Fy(mapD));

% set up Neumann boundary conditions
qN = zeros(Nfp*Nfaces, K);
%qN(mapN) = nx(mapN).*(pi*cos(pi*Fx(mapN)).*sin(pi*Fy(mapN))) + ...
%           ny(mapN).*(pi*sin(pi*Fx(mapN)).*cos(pi*Fy(mapN))) ;

% evaluate boundary condition contribution to rhs
Aqbc = PoissonIPDGbc2D(uD, qN);

% set up right hand side forcing
%rhs = -2*(pi^2)*sin(pi*x).*sin(pi*y);
[~, rhs] = solution1(x,y);
rhs = MassMatrix*(J.*rhs) + Aqbc;

% % solve system
% u = A\rhs(:);
% u = reshape(u, Np, K);
% SurfDG2D(u)

% % For plotting use
%[TRI,xout,yout,uout,interp] = PlotField2D(2*N,x,y,u);
%clf;
%PlotDomain2D();
%hold on;
%PlotContour2D(TRI,xout,yout,uout,linspace(0,-12000,12));
%hold off;

%spy(A);

% sep_tree = GenNestedDissection2D(nmax);
% sep_tree = GenNestedDissectionBoxes2D(nmax);
%elim_tree = GenElimTreeOld2D(nmax);
elim_tree = GenElimTree2D(nmax, 'PartitionMode', 'alternating');
% PlotSparsityReorder(A, sep_tree);
%[perm, iperm] = ElimTree2Perm(elim_tree, 'postorder');
%figure
%spy(A(perm, perm))
%VisualizeTree2D(elim_tree)
elim_tree = ElimTree2SequentialTree(elim_tree);


% Save matrix, RHS and sep_tree in binary file
% Matlab stores matrices in CSC format
% If we want CSR format we can transpose before saving
% A = A'; transposed = true;
transposed = false;
b = rhs(:);
save(filename,'A','b','elim_tree','transposed','-v7.3','-nocompression');

return
end

% define analytical solutions
function [u, f] = solution1(x,y)
  u = (1-x.^2)/2;
  for k=1:2:23
    u = u - 16/pi^3 * sin(k*pi*(1+x)/2)/(k^3*sinh(k*pi)).*( sinh(k*pi*(1+y)/2) + sinh(k*pi*(1-y)/2));
  end
  f = ones(size(x));
end
