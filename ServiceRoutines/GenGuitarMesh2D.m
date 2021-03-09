function [VX, VY, K, EToV, BC] = GenGuitarMesh2D(h0, varargin)

% function [VX, VY, K, EToV, BCType] = GenSquareQuadMesh2D(h0)
% Purpose  : Generate 2D guitar mesh;   
% written by Boris Bonev

bl0 = 35/h0;

BCType = 'Dirichlet'; % Use Dirichlet by default
if nargin > 1
    if strcmp(varargin{1}, 'Dirichlet') || strcmp(varargin{1}, 'Neumann')...
        || strcmp(varargin{1}, 'mixed')
      BCType = varargin{1};
    else
      error('GenGuitarMesh2D: Unknown boundary type specified.')
    end
end

fd=@(p) dguitar(p);
fh = @huniform;
Bbox = [-8,-19.675;8,20];
pfix = [1,0;-1,0;0,-19.675;-0.6,-19.675;0.6,-19.675];

% Call distmesh
%[Vert,EToV]=distmesh2d(fd,fh,bl0,Bbox,pfix);
[Vert,EToV]=mydistmesh2d(fd,fh,bl0,Bbox,pfix,500);
VX = Vert(:,1); VY = Vert(:,2);
Nv = length(VX); K  = size(EToV,1);


% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);

% Boundary conditions
[EToE,EToF]= tiConnect2D(EToV);

% initialize BC flags
BC = zeros(K,3);

% loop over elements
for i=1:K
    
  % loop over faces
  for j=1:3
        
    % if no neighbor
    if (EToE(i,j) == i)
      if strcmp(BCType, 'Dirichlet')
        BC(i,j) = 6;
      elseif strcmp(BCType, 'Neumann')
        BC(i,j) = 7;
      elseif strcmp(BCType, 'mixed')
        if ( -1.9500 >= VY(EToV(i,j)) ) && ( VY(EToV(i,j)) >= -10.75 )...
            && ( -4 <= VX(EToV(i,j)) ) && ( VX(EToV(i,j)) <= 4 )
          BC(i,j) = 7;
        else
          BC(i,j) = 6;
        end
      end
    end
  end
end

% flatten the coordinate vectors
VX = VX(:)';  VY = VY(:)';

return
end

% x = linspace(-5,19.675,200);
% y = linspace(-8,8,100);
% [X,Y] = meshgrid(x,y);
% Z = Y.^2 - gibson_j40(X).^2;
% %Z(X < 0) = -X(X < 0);
% surf(X,Y,Z)

function dist=dguitar(p)
   p = protate(p,1.5*pi);
   dist = ddiff(dunion(dguitar_body(p), drectangle(p,-20,1,-1,1)), dcircle(p,6.25,0,2.4));
end

function dist=dguitar_body(pc)

x = 19.675*0.5*(JacobiGL(0,0,500)+1);
%x = linspace(0,19.675,1000);
y = gibson_j40(x);
z = - y;
x = [x, x(end:-1:1)];
y = [y, z(end:-1:1)];

dist = zeros(size(pc,1),1);
for i=1:size(pc,1)
  out = abs(pc(i,2)) > abs(gibson_j40(pc(i,1))) || pc(i,1) >= 19.675 || pc(i,1) < 0;
  dd = sqrt((pc(i,1)-x).^2+(pc(i,2)-y).^2);
  if out
    dist(i) = min(dd,[],'all');
  else
    dist(i) = -1*min(dd,[],'all');
  end
end

% dist = pc(:,2).^2 - gibson_j40(pc(:,1)).^2;
% sel = pc(:,1) <= 0;
% dist(sel) = inf;%pc(sel,1).^2 + pc(sel,2).^2 ;
% sel = pc(:,1) >= 19.675;
% dist(sel) = inf;%(pc(sel,1)-19.675).^2 + pc(sel,2).^2 ;
end

function y = gibson_j40(x)

a = 1.185387;
b = 0.3756804;
c = 2.940726;
d = 8.62169;
e = 0.08442539;
f = -0.5786059;
g = 5.001006;
h = 0.1806989;
i = 1.527422;
j = 0.05155868;
k = -0.9508379;
l = -0.1762013;
m = 3.092091;
n = 0.5851531;
o = -0.1935072;
p = 0.01593922;
q = -0.000402921;
r = 21.81638;
s = 34.6247;
bl = 19.675;

y = (a*sin(b*x + c) + d*sin(e*x + f) + g*sin(h*x + i) + j*sin(k*x + l)).*...
  (m + n*x + o*x.^2 + p*x.^3 + q*x.^4).*atan(r*x).*atan(s*(19.675 - x));
end
