l2err = [];
pvals = [];
hvals = [];

for h = 2.^(2:5)
for N = 1:3
  % Driver script for solving the 2D Poisson equation
  Globals2D;

  % generate the mesh
  [VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h, h, 'Dirichlet');
  %[VX, VY, K, EToV] = GenCircleMesh2D(1/h);

  % Initialize DG solver and construct grid and metric
  StartUp2D;

  % Compute translation maps that map DG data onto CG datastructures
  BuildGlobalMaps2D;
  
  % Construct analytical solution
  [uA, f] = solution1(x,y);

  % Initialize solver and construct grid and metric
  [A,M] = PoissonIPDG2D();
  
  % set up boundary conditions
  uD = zeros(Nfp*Nfaces, K);
  qN = zeros(Nfp*Nfaces, K);
  
  % compute rhs
  Aqbc = PoissonIPDGbc2D(uD, qN);
  rhs = MassMatrix*(J.*f) + Aqbc;

  % solve the numerical proble
  u = A\rhs(:);
  uA = uA(:);
  
  % compute L2 error
  e = u - uA;
  % comparison
  close all
  figure
  SurfDG2D(reshape(u,Np,K));
  figure
  SurfDG2D(reshape(uA,Np,K));
  figure
  SurfDG2D(reshape(e,Np,K));
  
  l2err = [l2err, sqrt( e'*M*e/(uA'*M*uA) )];
  pvals = [pvals, N];
  hvals = [hvals, h];
end
end

figure
for N=unique(pvals)
  selector = (pvals == N);
  loglog(1./hvals(selector), l2err(selector), '-o', 'DisplayName', sprintf('p = %d', N))
  hold on
  loglog(1./hvals(selector), hvals(selector).^-(N+1), 'DisplayName', sprintf('O(h^%d)', N+1))
  hold on
end
title('Convergence behavior for different polynomial orders')
xlabel('mesh width 1/h');
ylabel('relative error');
set(gca, 'xdir', 'reverse');
legend('Location','northwest');
hold off


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