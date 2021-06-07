function ContfdDG2D(u, vx, vy, varargin)

% function ContfdDG2D(u, vx, vy)
% Purpose: Plot deformedDG solution. elements are deformed according to the
% deformation field

  Globals2D;
  
  eps = 1E-5;
  if nargin > 3
    alpha = varargin{1};
  else
    alpha = 1E-3;
  end

  xb = ones(size(x,1),1)*(sum(x,1)/size(x,1)); yb = ones(size(y,1),1)*sum(y,1)/size(y,1);
  % apply minimal separation
  xDG = xb + (1-eps)*(x-xb); yDG = yb + (1-eps)*(y-yb);
  % compute triangulation
  dt = delaunayTriangulation(xDG(:),yDG(:));
  tri = dt.ConnectivityList;
  % apply deformation
  xDG = xDG + vx*alpha;
  yDG = yDG + vy*alpha;
  %figure
  %trimesh(tri,xDG,yDG,u(:))
  %h=trisurf(tri,xDG,yDG,u(:));
  %set(h,'LineWidth',0.01)
  [CS,h]=tricontf(xDG(:),yDG(:),tri,u(:), linspace(min(u(:)), max(u(:)), 20));
  %set(h,'edgecolor','none');
  set(h,'LineWidth',0.1);
  %set(h,'LineColor','none')
  colormap(parula(20))
  colorbar
  %shading interp
return;