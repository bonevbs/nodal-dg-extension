function SurfDG2D(u, varargin)

% function SurfDG2D(u)
% Purpose: Plot DG solution. elements are scaled back by epsilon to allow
%          the triangulation to work.

cont = false; % contour lines by default
sym_scale = false;
if nargin > 1
    if varargin{1}
      cont = true;
    else
      cont = false;
    end
end

Globals2D;

Nout = N;

%epsilon = 0;
epsilon = 1E-5;

xb = ones(size(x,1),1)*(sum(x,1)/size(x,1)); yb = ones(size(y,1),1)*sum(y,1)/size(y,1);
xDG = xb + (1-epsilon)*(x-xb); yDG = yb + (1-epsilon)*(y-yb);

[tri, xout, yout, uout, interp] = triangulate(Nout, x, y, u);

%dt = delaunayTriangulation(xDG(:),yDG(:));
%tri = dt.ConnectivityList;
trisurf(tri, xout(:), yout(:), uout(:));
shading interp
if sym_scale
  sc = max(abs(uout(:)));
  caxis([-sc sc])
end
colorbar('southoutside')
view([0 90])
hold on
%plot_domain(u);
hold on
if cont; contour_lines(tri, xout(:), yout(:), uout(:), linspace(-sc, sc, 11)); end
hold off


%figure
%trimesh(tri,xDG,yDG,u(:))
%h=trisurf(tri,xDG,yDG,u(:));
%set(h,'LineWidth',0.01)
%[CS,h]=tricontf(xout(:), yout(:), TRI, uout(:), linspace(min(uout(:)), max(uout(:)), 20));
%set(h,'edgecolor','none');
%set(h,'LineWidth',0.1);
%set(h,'LineColor','none')
%colormap(parula(20))
%colorbar
%shading interp
return;
end

function [TRI, xout, yout, uout, interp] = triangulate(Nout, xin, yin, uin)

Globals2D;

% build equally spaced grid on reference triangle
Npout = (Nout+1)*(Nout+2)/2;
rout = zeros(Npout,1); sout = zeros(Npout,1);
sk = 1;
for n=1:Nout+1
  for m=1:Nout+2-n
    rout(sk) = -1 + 2*(m-1)/Nout;
    sout(sk) = -1 + 2*(n-1)/Nout;
    counter(n,m) = sk; sk = sk+1;
  end
end

% build matrix to interpolate field data to equally spaced nodes
interp = InterpMatrix2D(rout, sout);

% build triangulation of equally spaced nodes on reference triangle
tri = [];
for n=1:Nout+1
  for m=1:Nout+1-n
    v1 = counter(n,m);   v2 = counter(n,m+1);
    v3 = counter(n+1,m); v4 = counter(n+1,m+1);
    if(v4)
      tri = [tri;[[v1 v2 v3];[v2 v4 v3]]];
    else
      tri = [tri;[[v1 v2 v3]]];
    end
  end
end

% build triangulation for all equally spaced nodes on all elements
TRI = [];
for k=1:K
  TRI = [TRI; tri+(k-1)*Npout];
end

% interpolate node coordinates and field to equally spaced nodes
xout = interp*xin; yout = interp*yin; uout = interp*uin;
end


function contour_lines(tri, x, y, u, levels)

% Purpose: generic routine to plot contours for triangulated data

draw_z = max(u, [], 'all');
  
Nlevels = length(levels);

hold on;

v1 = tri(:,1); v2 = tri(:,2); v3 = tri(:,3);
u1 = u(v1);    u2 = u(v2);    u3 = u(v3);
x1 = x(v1);    x2 = x(v2);    x3 = x(v3);
y1 = y(v1);    y2 = y(v2);    y3 = y(v3);

allx = []; ally = [];
for n=1:Nlevels
  lev = levels(n);

  flag1 = max(u1,u2)>=lev & min(u1,u2)<=lev;   % edge 1
  flag2 = max(u2,u3)>=lev & min(u2,u3)<=lev;   % edge 2
  flag3 = max(u1,u3)>=lev & min(u1,u3)<=lev;   % edge 3
  
  c1 = (lev-u1)./(u2-u1);  xc1 = (1-c1).*x1 + c1.*x2; yc1 = (1-c1).*y1 + c1.*y2;
  c2 = (lev-u2)./(u3-u2);  xc2 = (1-c2).*x2 + c2.*x3; yc2 = (1-c2).*y2 + c2.*y3;
  c3 = (lev-u1)./(u3-u1);  xc3 = (1-c3).*x1 + c3.*x3; yc3 = (1-c3).*y1 + c3.*y3;
  
  ids = find(flag1+flag2==2);
  allx = [allx, [xc1(ids),xc2(ids)]'];
  ally = [ally, [yc1(ids),yc2(ids)]'];

  ids = find(flag2+flag3==2);
  allx = [allx, [xc2(ids),xc3(ids)]'];
  ally = [ally, [yc2(ids),yc3(ids)]'];

  ids = find(flag1+flag3==2);
  allx = [allx, [xc1(ids),xc3(ids)]'];
  ally = [ally, [yc1(ids),yc3(ids)]'];
end
hold off

ha = line(allx, ally, draw_z*ones(size(allx)));
set(ha, 'Color', 'black')

hold off;
return;
end

function plot_domain(u)

% Purpose: Show domain boundary

draw_z = max(u, [], 'all');

Globals2D;

hold on;

for k=1:K
  for f=1:Nfaces
    bc = BCType(k,f);
    if bc ~= 0
      ids = (k-1)*Nfp*Nfaces+(f-1)*Nfp+(1:Nfp);
      ha = line(Fx(ids), Fy(ids), draw_z*ones(size(Fx(ids))));
      set(ha, 'Color', 'black')
    end
    ids = (k-1)*Nfp*Nfaces+(f-1)*Nfp+(1:Nfp);
  end
end

hold off
return;
end
