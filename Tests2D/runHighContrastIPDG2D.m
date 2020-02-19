clear all

% Driver script for solving the 2D Poisson equation
Globals2D;
GlobalsCG;

% Polynomial order used for approximation
N=1;

% generate the mesh
[VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(10, 10);
%[VX, VY, K, EToV] = GenCircleMesh2D(1/4);

% Initialize solver and construct grid and metric
StartUp2D;

% set up boundary conditions
BuildBCMaps2D;

% set up Dirichlet boundary conditions
uD = zeros(Nfp*Nfaces, K);
[uD(mapD),~,aD] = solution2(Fx(mapD),Fy(mapD));

% set up Neumann boundary conditions
qN = zeros(Nfp*Nfaces, K);
aN = c*ones(Nfp*Nfaces, K);
%qN(mapN) = nx(mapN).*(pi*cos(pi*Fx(mapN)).*sin(pi*Fy(mapN))) + ...
%           ny(mapN).*(pi*sin(pi*Fx(mapN)).*cos(pi*Fy(mapN))) ;

% Compute the stiffness and mass matrices
%[A,M] = PoissonIPDG2D();
[A,M] = HighContrastIPDG2D(a);

% evaluate boundary condition contribution to rhs
Aqbc = HighContrastIPDGbc2D(uD, qN, a, aN);

% set up right hand side forcing
[uExact,rhs] = solution2(x,y);
rhs = c*MassMatrix*(J.*rhs) + Aqbc;

% solve system
u = A\rhs(:);
u = reshape(u, Np, K);

e = u(:) - uExact(:);
fprintf("Relative L2 error: %e\n", sqrt( e'*M*e/(uExact(:)'*M*uExact(:)) ));

% do some plotting
[X,Y] = meshgrid(linspace(-1,1,20));
Z = griddata(x,y,u,X,Y,'cubic');
surf(X,Y,Z)
%dt = delaunayTriangulation(x,y);
%tri = dt.ConnectivityList ;
%figure
%trisurf(tri,x,y,u)
%figure
%scatter3(xCG, yCG, uCG)
%uDG = zeros(K,Np);
%uDG = uCG(gmap)';
%clf
%PlotMesh2D;
%hold on
%PlotAdaptiveContour2D(uDG,linspace(0,1,12),1e-1);
%hold off

% define analytical solutions
function [u, f] = solution1(x,y)
  u = (1-x.^2)/2;
  for k=1:2:23
    u = u - 16/pi^3 * sin(k*pi*(1+x)/2)/(k^3*sinh(k*pi)).*( sinh(k*pi*(1+y)/2) + sinh(k*pi*(1-y)/2));
  end
  f = ones(size(x));
end

function [u, f] = solution2(x,y)
  u = x.^2 + y.^2;
  f = -4*ones(size(x));
end
