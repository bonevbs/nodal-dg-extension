l2err = [];
pvals = [];
hvals = [];

for h = 2.^(2:7)
for N = 1:4
  % Driver script for solving the 2D High Contrast Poisson equation
  Globals2D;

  % generate the mesh
  [VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h, h);

  % Initialize solver and construct grid and metric
  StartUp2D;

  % set up boundary conditions
  BuildBCMaps2D;

  % Checkerboard pattern for solution 3
  a1 = 4; a2 = 1;
  EP = (1:K)';
  E1 = find( ( VX(EToV(EP, 1)) + VX(EToV(EP, 2)) + VX(EToV(EP, 3)) )/3 .* ( VY(EToV(EP, 1)) + VY(EToV(EP, 2)) + VY(EToV(EP, 3)) )/3 <= 0)';
  E2 = find( ( VX(EToV(EP, 1)) + VX(EToV(EP, 2)) + VX(EToV(EP, 3)) )/3 .* ( VY(EToV(EP, 1)) + VY(EToV(EP, 2)) + VY(EToV(EP, 3)) )/3 > 0)';
  a = ones(K,1);
  a(E1) = a1*a(E1);
  a(E2) = a2*a(E2);
  aElements = ones(size(x,1),1)*a';
  aFaces = aElements(Fmask(:), :);

  
  % Compute the stiffness and mass matrices
  [A,M] = HighContrastIPDG2D(a);

  % set up Dirichlet boundary conditions
  uD = zeros(Nfp*Nfaces, K);
  [uD(mapD),~] = solution3(Fx(mapD),Fy(mapD),aFaces(mapD));

  % set up Neumann boundary conditions
  qN = zeros(Nfp*Nfaces, K);
  aN = ones(Nfp*Nfaces, K);

  % evaluate boundary condition contribution to rhs
  Aqbc = HighContrastIPDGbc2D(uD, qN, a, aN);

  % set up right hand side forcing and adjust according to the jump in a
  [uExact,f] = solution3(x,y,aElements);
  rhs = MassMatrix*(J.*f) + Aqbc;

  % solve system
  u = A\rhs(:);
  u = reshape(u, Np, K);
  
  % compute L2 error
  e = u(:) - uExact(:);
  l2err = [l2err, sqrt( e'*M*e/(uExact(:)'*M*uExact(:)) )];
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
function [u, f] = solution3(x,y,a)
  u = sin(pi*x).*sin(pi*y)./a;
  f = 2*pi^2*sin(pi*x).*sin(pi*y);
end