clear all

% Driver script for solving the 2D Poisson equation
Globals2D;
GlobalsCG;

% Polynomial order used for approximation
N=25;

% generate the mesh
[VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(1, 1);
%[VX, VY, K, EToV] = GenCircleMesh2D(1/4);

% Initialize DG solver and construct grid and metric
StartUp2D;

% Compute translation maps that map DG data onto CG datastructures
BuildGlobalMaps2D;

% Compute the stiffness and mass matrices
[A,M] = PoissonCG2D();

% compute grid points in CG enumeration
xCG = zeros(nTotal,1); yCG = zeros(nTotal,1);
xCG(gmap') = x; yCG(gmap') = y;

% compute exact solution for boundary conditions
[uExact, fExact] = solution1(xCG,yCG);

% solve the numerical proble
u = A(1:nFree, 1:nFree) \ ( M(1:nFree, 1:nTotal)*fExact - A(1:nFree, nFree+1:nTotal)*uExact(nFree+1:nTotal) );
uCG = [u; uExact(nFree+1:nTotal)];

eCG = uCG - uExact;
fprintf("Relative L2 error: %e\n", sqrt( eCG'*M*eCG/(uExact'*M*uExact) ));

% do some plotting
SurfCG2D(uCG)

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
