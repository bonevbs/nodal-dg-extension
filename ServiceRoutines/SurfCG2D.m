function SurfCG2D(u)

% function SurfCG2D(u)
% Purpose: Plot CG solution. gmap is used to infer mesh structure.

  Globals2D;
  GlobalsCG;

  % compute grid points in CG enumeration
  xCG = zeros(nTotal,1); yCG = zeros(nTotal,1);
  xCG(gmap') = x; yCG(gmap') = y;

  dt = delaunayTriangulation(xCG(:),yCG(:));
  tri = dt.ConnectivityList;
  figure
  trimesh(tri,xCG,yCG,u(:))
return;