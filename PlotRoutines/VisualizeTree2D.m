function VisualizeTree2D(elim_tree)

Globals2D;

eps = 1E-5;

v = zeros(Np*K,1);
for i=1:size(elim_tree,2)
  v(elim_tree{1,i}) = i;
end
v=reshape(v,Np,K);


xb = ones(size(x,1),1)*(sum(x,1)/size(x,1)); yb = ones(size(y,1),1)*sum(y,1)/size(y,1);
xDG = xb + (1-eps)*(x-xb); yDG = yb + (1-eps)*(y-yb);

%dt = delaunayTriangulation(xDG(:),yDG(:));
%tri = dt.ConnectivityList;
[tri, xout, yout, vout, interp] = triangulate(N, x, y, v);

trisurf(tri,xout(:),yout(:),vout(:))
colormap( parula(size(elim_tree,2)) )
shading flat
view(0,90)
hold on
%plot_domain();
hold off

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

function plot_domain()

% Purpose: Show domain boundary

Globals2D;

hold on;

for k=1:K
  for f=1:Nfaces
    bc = BCType(k,f);
    if bc ~= 0
      ids = (k-1)*Nfp*Nfaces+(f-1)*Nfp+(1:Nfp);
      ha = line(Fx(ids), Fy(ids), 1000*ones(size(Fx(ids))));
      set(ha, 'Color', 'black')
    end
    ids = (k-1)*Nfp*Nfaces+(f-1)*Nfp+(1:Nfp);
  end
end

hold off
return;
end