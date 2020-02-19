l2err = [];
pvals = [];
hvals = [];

for h = 2.^(2:4)
for N = 1:3
  % Driver script for solving the 2D Poisson equation
  Globals2D;
  GlobalsCG;

  % generate the mesh
  [VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h, h);
  %[VX, VY, K, EToV] = GenCircleMesh2D(1/h);

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
  
  % compute L2 error
  eCG = uCG - uExact;
  l2err = [l2err, sqrt( eCG'*M*eCG/(uExact'*M*uExact) )];
  pvals = [pvals, N];
  hvals = [hvals, h];
end
end

figure
for N=unique(pvals)
  selector = (pvals == N);
  loglog(1./hvals(selector), l2err(selector), 'DisplayName', sprintf('p = %d', N))
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
