clear all

l2err = [];
pvals = [];
hvals = [];

icase = 'forced 1';

for h = 2.^(2:6)
  for p = 1:4
    clear global
    
    % Driver script for solving the 2D Poisson equation
    Globals2D;
    N = p;
    
    % generate the mesh
    if strcmp(icase, 'forced 1')
      [VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h, h, 'Dirichlet');
    elseif strcmp(icase, 'forced 2')
      [VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h, h, 'channel');
    else
      error("Unknown geometry.")
      %[VX, VY, K, EToV, BCType] = GenCircleMesh2D(h0);
    end
    
    % Initialize DG solver and construct grid and metric
    StartUp2D;
    
    % Build boundary maps
    BuildBCMaps2D;
    
    % set up Dirichlet boundary conditions
    uxD = zeros(Nfp*Nfaces, K);
    uyD = zeros(Nfp*Nfaces, K);
    % set up Neumann boundary conditions
    qxN = zeros(Nfp*Nfaces, K);
    qyN = zeros(Nfp*Nfaces, K);
    
    switch icase
      case 'forced 1'
        sel1 = ( VX(EToV(:, 1)) + VX(EToV(:, 2)) + VX(EToV(:, 3)) )/3 > 0.5;
        [uxA, uyA, fx, fy, mu, lambda] = solution1(x,y,sel1);
        uxD(mapD) = uxA(vmapD);
        uyD(mapD) = uyA(vmapD);
      case 'forced 2'
        [uxA, uyA, qxA, qyA, fx, fy, mu, lambda] = solution2(x,y);
        uxD(mapD) = uxA(vmapD);
        uyD(mapD) = uyA(vmapD);
        qxN(mapN) = qxA(vmapN);
        qyN(mapN) = qyA(vmapN);
    end
    
    % Initialize solver and construct grid and metric
    [Axx,Axy,Ayx,Ayy,M] = ElasticityIPDG2D(mu, lambda);
    A = [Axx, Axy; Ayx, Ayy];
    
    % evaluate boundary condition contribution to rhs
    [bcx, bcy] = ElasticityIPDGbc2D(uxD,uyD,qxN,qyN);
    
    % set up right hand side forcing and adjust according to the jump in a
    fx = - MassMatrix*(J.*fx); fy = - MassMatrix*(J.*fy);
    rhs = [fx(:); fy(:)] + [bcx(:); bcy(:)];
    
    % create big Mass matrix
    M = blkdiag(M,M);
    
    % solve the numerical proble
    u = A\rhs(:);
    uA = [uxA(:); uyA(:)];
    
    % compute L2 error
    e = u - uA;
    % comparison
    close all
    %figure
    %SurfDG2D(reshape(lambda,Np,K));
    %figure
    %SurfDG2D(reshape(mu,Np,K));
    figure
    SurfDG2D(reshape(u(1:Np*K),Np,K));
    figure
    SurfDG2D(reshape(uxA,Np,K));
    figure
    SurfDG2D(reshape(e(1:Np*K),Np,K));
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
function [ux, uy, fx, fy, mu, lambda] = solution1(x,y,sel1)

c = 0;

muc0 = 1; mux0 = 0; muy0 = 0;
muc1 = 10000; mux1 = 0; muy1 = 0;
lambdac0 = 1; lambdax0 = 0; lambday0 = 0;
lambdac1 = 10000; lambdax1 = 0; lambday1 = 0;
ux = sin(pi*x).*sin(pi*y) + c;
uy = sin(pi*x).*sin(pi*y) + c;
mu = muc0 + mux0*x + muy0*y;
lambda = lambdac0 + lambdax0*x + lambday0*y;
mu(:,sel1) = muc1 + mux1*x(:,sel1) + muy1*y(:,sel1);
mux = mux0 * ones(size(x));
mux(:,sel1) = mux1;
muy = muy0 * ones(size(x));
muy(:,sel1) = muy1;
lambda(:,sel1) = lambdac1 + lambdax1*x(:,sel1) + muy1*y(:,sel1);
lambdax = lambdax0 * ones(size(x));
lambdax(:,sel1) = lambdax1;
lambday = lambday0 * ones(size(x));
lambday(:,sel1) = lambday1;

fx = pi*lambdax.*cos(pi*x).*sin(pi*y) + pi*lambdax.*sin(pi*x).*cos(pi*y)...
  - pi^2*lambda.*sin(pi*x).*sin(pi*y) + pi^2*lambda.*cos(pi*x).*cos(pi*y)...
  + 2*pi*mux.*cos(pi*x).*sin(pi*y) - 2*pi^2*mu.*sin(pi*x).*sin(pi*y)...
  + pi*muy.*cos(pi*x).*sin(pi*y) + pi*muy.*sin(pi*x).*cos(pi*y)...
  - pi^2*mu.*sin(pi*x).*sin(pi*y) + pi^2*mu.*cos(pi*x).*cos(pi*y);
fy = pi*lambday.*cos(pi*x).*sin(pi*y) + pi*lambday.*sin(pi*x).*cos(pi*y)...
  - pi^2*lambda.*sin(pi*x).*sin(pi*y) + pi^2*lambda.*cos(pi*x).*cos(pi*y)...
  + 2*pi*muy.*sin(pi*x).*cos(pi*y) - 2*pi^2*mu.*sin(pi*x).*sin(pi*y)...
  + pi*mux.*cos(pi*x).*sin(pi*y) + pi*mux.*sin(pi*x).*cos(pi*y)...
  - pi^2*mu.*sin(pi*x).*sin(pi*y) + pi^2*mu.*cos(pi*x).*cos(pi*y);
end

function [uxA, uyA, qxA, qyA, fx, fy, mu, lambda] = solution2(x,y)
uxA = sin(pi*x).*cosh(pi*y);
uyA = -cos(pi*x).*sinh(pi*y);
lambda = 1.0*ones(size(x));
mu = 1.0*ones(size(x));
fx = zeros(size(x));
fy = zeros(size(x));
qxA = 2*pi*sign(y).*mu.*sin(pi*x).*sinh(pi*y);
qyA = -2*pi*sign(y).*mu.*cos(pi*x).*cosh(pi*y);
end
