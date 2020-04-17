%====================================================================================
% -- IN --
% numx,y - number of subdivisions along x- and y-axis
%
% -- OUT --
% Nv     - number of nodes
% VX     - x-coordinates of nodes
% VY     - y-coordinates of nodes
% K      - number of elements
% EToV   - elemet-to-vertices connectivities
% BC    - boundary conditions flags
%====================================================================================


function [Nv, VX, VY, K, EToV, BC] = GenStructuredSquareQuadMesh2D(numx,numy)

% numx,numy - number of subdivisions along x- and y-axis

% Nodes enumeration
%
%   13 --14 --15 --16
%    |    |    |    |
%    9 --10 --11 --12
%    |    |    |    |
%    5 -- 6 -- 7 -- 8
%    |    |    |    |
%    1 -- 2 -- 3 -- 4

% x- and y-coordinate of nodes
x = (0:1/(numx):1);
y = (0:1/(numy):1);

% Matlab's meshgrid is used to create 2D grid from specified divisions above
[X,Y] = meshgrid(x,y);

% number of vertices
Nv = length(x)*length(y);

% reshape into 1-D arrays
VX = reshape(X',Nv,1);
VY = reshape(Y',Nv,1);

% Build up quad-to-vertices connectivities
element = zeros(numx*numy,3);

k=0;
for j=1:numy
    for i=1:numx
        k = k+1;
        
        element(k,1) = (j-1)*(numx+1) + i;
        element(k,2) = (j-1)*(numx+1) + i+1;
        element(k,3) =  j   *(numx+1) + i+1;
        element(k,4) =  j   *(numx+1) + i;
        
    end
end

% generate triangle-to-vertices connectivities
EToV = triangulate(VX,VY,element);

VX = VX';
VY = VY';

% number of triangular elements
K = size(EToV,1);



% -- BOUNDARY CONDITIONS --------------------------------------------------
[EToE,EToF]= tiConnect2D(EToV);

% Dirichlet BC flag
Dirichlet = 6;

% initialize BC flags
BC = zeros(K,3);

% loop over elements
for i=1:K
    
    % loop over faces
    for j=1:3
        
        % if no neighbor
        if (EToE(i,j) == i)
            BC(i,j) = Dirichlet;
        end
    end
end





% % 
% % Z= zeros(length(y)-2,length(x)-2);
% % 
% nn = EToV; % nodes in a long list
% 
% xx = VX(nn); yy = VY(nn); % coordinates corresponding to nodes
% 
% xplot = reshape(xx,size(EToV));
% 
% yplot = reshape(yy,size(EToV));figure(1);
% 
% clf;
% 
% fill(xplot',yplot','w');
% 
% title(['Triangular Mesh ', num2str(length(nn)),' elements'])

end


% service function
% ================

% -- PURPOSE --
% Given quad-to-vertices connectivities "nodes", generate 
% triangle-to-vertices connectivities "tri" by shortest diagonal
% refinements. This requires knowledge of the x-coordinates "xnod" and
% y-coordinates "ynod" of the vertices.

function tri = triangulate(xnod,ynod,nodes)

    %
    % triangulate the quadrilateral mesh
    %
    
    nele = size(nodes,1);
    tri = zeros(3,2*nele)';
    
    ii1 = [2 3 1];
    jj1 = [4 1 3];
    ii2 = [1 2 4];
    jj2 = [2 3 4];
    
    % number of trinagles created so far
    nrtri = 0;
    for iel = 1:nele
        
        %  4 ----- 3        4 ----- 3
        %  | \     |        |     / |
        %  |   \   |   OR   |   /   |
        %  |     \ |        | /     |
        %  1 ----- 2        1 ----- 2
        %     d2               d1
        
        iv = nodes(iel,:);
        d1 = norm([xnod(iv(1))-xnod(iv(3)) ; ynod(iv(1))-ynod(iv(3))]);
        d2 = norm([xnod(iv(2))-xnod(iv(4)) ; ynod(iv(2))-ynod(iv(4))]);
        
        if d1 <= d2
            nrtri = nrtri+1;
            tri(nrtri,:) = iv(ii1);
            nrtri = nrtri+1;
            tri(nrtri,:) = iv(jj1);
        else
            nrtri = nrtri+1;
            tri(nrtri,:) = iv(ii2);
            nrtri = nrtri+1;
            tri(nrtri,:) = iv(jj2);
        end
    end
end