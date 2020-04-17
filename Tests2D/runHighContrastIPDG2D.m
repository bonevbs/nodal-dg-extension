clear all

% Driver script for solving the 2D Poisson equation
Globals2D;

% Polynomial order used for approximation
N=2;

% generate the mesh
[VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(20, 20);
%[VX, VY, K, EToV] = GenCircleMesh2D(1/4);

% Initialize solver and construct grid and metric
StartUp2D;

% set up boundary conditions
BuildBCMaps2D;

% Checkerboard pattern for solution 3
EP = (1:K)';
E1 = find( ( VX(EToV(EP, 1)) + VX(EToV(EP, 2)) + VX(EToV(EP, 3)) )/3 .* ( VY(EToV(EP, 1)) + VY(EToV(EP, 2)) + VY(EToV(EP, 3)) )/3 <= 0)';
E2 = find( ( VX(EToV(EP, 1)) + VX(EToV(EP, 2)) + VX(EToV(EP, 3)) )/3 .* ( VY(EToV(EP, 1)) + VY(EToV(EP, 2)) + VY(EToV(EP, 3)) )/3 > 0)';

a1 = 1; a2 = 20;
a = ones(K,1);
a(E1) = a1*a(E1);
a(E2) = a2*a(E2);
Ea = ones([size(x,1)],1)*a';
Fa = Ea(Fmask(:), :);


% Compute the stiffness and mass matrices
%[A,M] = PoissonIPDG2D();
[A,M] = HighContrastIPDG2D(a);

% set up Dirichlet boundary conditions
uD = zeros(Nfp*Nfaces, K);
%[uD(mapD),~] = solution1(Fx(mapD),Fy(mapD));
[uD(mapD),~] = solution3(Fx(mapD),Fy(mapD),Fa(mapD));

% set up Neumann boundary conditions
qN = zeros(Nfp*Nfaces, K);
%aN = ones(Nfp*Nfaces, K);

% evaluate boundary condition contribution to rhs
Aqbc = HighContrastIPDGbc2D(uD, qN, a);

% set up right hand side forcing and adjust according to the jump in a
%[uExact,rhs] = solution1(x,y);
[uExact,f] = solution3(x,y,Ea);
rhs = MassMatrix*(J.*f) + Aqbc;

% solve system
u = A\rhs(:);
u = reshape(u, Np, K);

e = u(:) - uExact(:);
fprintf("Relative L2 error: %e\n", sqrt( e'*M*e/(uExact(:)'*M*uExact(:)) ));

% do some plotting
SurfDG2D(u);

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

function [u, f] = solution3(x,y,a)
  u = sin(pi*x).*sin(pi*y)./a;
  f = 2*pi^2*sin(pi*x).*sin(pi*y);
end
