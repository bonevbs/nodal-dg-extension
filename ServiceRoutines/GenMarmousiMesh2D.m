function [VX, VY, K, EToV, BC, Rho, Mu, Lambda] = GenMarmousiMesh2D(h0, varargin)

BCType = 'Dirichlet'; % Use Dirichlet by default
if nargin > 1
    if strcmp(varargin{1}, 'Dirichlet') || strcmp(varargin{1}, 'Neumann')...
        || strcmp(varargin{1}, 'corner') || strcmp(varargin{1}, 'channel')...
        || strcmp(varargin{1}, 'bottom') || strcmp(varargin{1}, 'Dirac')
      BCType = varargin{1};
    else
      error('GenMarmousiMesh2D: Unknown boundary type specified.')
    end
end


% load Marmousi data
load('/Users/boris/Code/Matlab/nodal-dg-extension/Grid2D/marmousi.mat');
% determine dimensions
xmin = min(xM, [], 'all');
xmax = max(xM, [], 'all');
ymin = min(yM, [], 'all');
ymax = max(yM, [], 'all');
ymax = -450; % cut off the water
h0 = (ymax - ymin)/h0;
Nx = round_even((xmax - xmin)/h0);
Ny = round_even((ymax - ymin)/h0);

% compute mu and lambda in the marmousi mesh
rho = 100^3 * 0.001 * rho;
mu = rho.*v_s.^2;
lambda = rho.*v_p.^2 - 2*mu;

% make an image
img = zeros(size(rho,1),size(rho,2),3);
img(:,:,1) = rho/max(rho,[],'all');
img(:,:,2) = v_p/max(v_p,[],'all');
img(:,:,3) = v_s/max(v_s,[],'all');
img = rgb2gray(img);
figure
imshow(img);

% segmentation using rho
thresh = multithresh(img,6);
seg_img = imquantize(img,thresh);
bnd_img = imgradient(seg_img);
figure
imshow(label2rgb(seg_img));


% generate mesh to sample points
mesh_img = zeros(size(seg_img));
[resy, resx] = size(mesh_img);
mesh_img(1:ceil(resy/Ny):end,:) = 1;
mesh_img(:,1:ceil(resx/Nx):end) = 1;
%mesh_img(1:ceil(1*resy/Ny):end,1:ceil(1*resx/Nx):end) = 1;
%figure
%imshow(mesh_img);

sel_img = bnd_img & mesh_img;
%figure
%imshow(sel_img);
xs = XM(sel_img);
ys = YM(sel_img);
pfix = [xs, ys];
pfix = remove_similar(pfix, 0.3*h0);
pfix = [pfix; xmin, ymin; xmin, ymax; xmax, ymin; xmax, ymax];
% create a uniform grid between 
xs = linspace(xmin, xmax, Nx);
ys = linspace(ymin, ymax, Ny);
[xs,ys] = meshgrid(xs, ys);
pnew = [xs(:), ys(:)];
% reject points based on image segmentation

p = [pfix; pnew];
p = remove_similar(p, 0.05*h0);
%scatter(pnew(:,1), pnew(:,2));

figure
scatter(p(:,1), p(:,2));
t = delaunay(p(:,1), p(:,2));
triplot(t,p(:,1), p(:,2));
[Vert,EToV]=fixmesh(p,t);

%mask = boundarymask(seg_img);
%RGB = label2rgb(seg_img);
%figure
%imshow(RGB)
%bnd_img = full(sprand(bnd_img) > 0.95);
%soft_edges = imgaussfilt(bnd_img,16);

% edge detection after quantization
%BW1 = edge(seg_img, 'log');
%imshow(BW1);

% method using contour lines
%figure
%[C,h] = imcontour(img,10);

% distmesh stuff
%fd = @(p) drectangle(p,xmin,xmax,ymin,ymax);
%[Vert,EToV] = mydistmesh2d(fd,@huniform,h0,[xmin,ymin;xmax,ymax],pfix,200);
VX = Vert(:,1); VY = Vert(:,2);
Nv = length(VX); K  = size(EToV,1);

% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);
VX = VX'; VY = VY';

% Boundary conditions
[EToE,EToF]= tiConnect2D(EToV);

% Dirichlet BC flag
Dirichlet = 6;
Neuman = 7;

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
      elseif strcmp(BCType, 'corner')
        if any(VX(EToV(i,1:3)) == xmin) || any(VY(EToV(i,1:3)) == ymin)
          BC(i,j) = 6;
        else
          BC(i,j) = 7;
        end
      elseif strcmp(BCType, 'channel')
        if ( VX(EToV(i,j)) == xmin && VX(EToV(i,mod(j,3)+1)) == xmin ) ||...
          ( VX(EToV(i,j)) == xmax && VX(EToV(i,mod(j,3)+1)) == xmax )
          BC(i,j) = 6;
        elseif ( VY(EToV(i,j)) == ymin && VY(EToV(i,mod(j,3)+1)) == ymin ) ||...
          ( VY(EToV(i,j)) == ymax && VY(EToV(i,mod(j,3)+1)) == ymax )
          BC(i,j) = 7;
        else
          error("Node should not be on the boundary.")
        end
      elseif strcmp(BCType, 'bottom')
        if ( VX(EToV(i,j)) <= xmin && VX(EToV(i,mod(j,3)+1)) <= xmin ) ||...
          ( VX(EToV(i,j)) >= xmax && VX(EToV(i,mod(j,3)+1)) >= xmax )
          BC(i,j) = 7;
        elseif ( VY(EToV(i,j)) >= ymax && VY(EToV(i,mod(j,3)+1)) >= ymax )
          BC(i,j) = 7;
        elseif ( VY(EToV(i,j)) <= ymin && VY(EToV(i,mod(j,3)+1)) <= ymin )
          BC(i,j) = 6;
        else
          error("Node should not be on the boundary.")
        end
      elseif strcmp(BCType, 'Dirac')
        if ( VX(EToV(i,j)) == (xmax-xmin)/2) || ( VY(EToV(i,j)) == ymax && VY(EToV(i,mod(j,3)+1)) == ymax )
          BC(i,j) = 6;
        else
          BC(i,j) = 7;
        end
      end
    end
  end
end

% compute barycentric coordinates for each cell
BX = ( VX(EToV(:,1)) + VX(EToV(:,2)) + VX(EToV(:,3)) )/3;
BY = ( VY(EToV(:,1)) + VY(EToV(:,2)) + VY(EToV(:,3)) )/3;

Rho = interp2(XM,YM,rho,BX,BY);
Mu = interp2(XM,YM,mu,BX,BY);
Lambda = interp2(XM,YM,lambda,BX,BY);

% flatten the coordinate vectors
VX = VX(:)';  VY = VY(:)';

return

end

function p = remove_similar(p, dist)
  n = size(p,1);
  xs = p(:,1)*ones(1,n);
  ys = p(:,2)*ones(1,n);
  
  rs = sqrt((xs - xs').^2 + (ys - ys').^2);
  sel = (rs < dist) - speye(size(rs));
  rem = [];
  for i=1:n
    if ~ismember(i, rem)
      torem = find(sel(i,:));
      torem = torem(torem > i);
      rem = union(rem, torem);
    end   
  end
  keep = setdiff(1:n, rem);
  p = p(keep,:);
end

function S = round_even(S)
% round to nearest even integer.
idx = mod(S,2)>1;
S = floor(S);
S(idx) = S(idx)+1;
end
