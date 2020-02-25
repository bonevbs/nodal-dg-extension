function SurfDG2D(u)

% function SurfDG2D(u)
% Purpose: Plot DG solution. elements are scaled back by epsilon to allow
%          the triangulation to work.

  Globals2D;
  
  eps = 1E-5;

  xb = ones(size(x,1),1)*(sum(x,1)/size(x,1)); yb = ones(size(y,1),1)*sum(y,1)/size(y,1);
  xDG = xb + (1-eps)*(x-xb); yDG = yb + (1-eps)*(y-yb);

  dt = delaunayTriangulation(xDG(:),yDG(:));
  tri = dt.ConnectivityList;
  figure
  trimesh(tri,xDG,yDG,u(:))
return;