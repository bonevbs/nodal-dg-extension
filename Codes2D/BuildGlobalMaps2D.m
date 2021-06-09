function [gmap, nTotal, nFree, nBoundary] = BuildGlobalMaps2D()

Globals2D;
GlobalsCG2D;

% this matrix stores nodes that are doubled in the form of a connectivity
% matrix
%doubleNodes = sparse(vmapM, vmapP, ones(K*Nfaces*Nfp,1));%+ sparse(vmapP, vmapM, ones(K*Nfaces*Nfp,1));
G = graph(vmapM, vmapP);

gmap = zeros(K*Np,1);
j = 1;
% unique counting to build a global map of degrees of freedom
for i=1:K*Np
  if gmap(i) == 0
    % go through all nodes corresponding to this node
    gmap(dfsearch(G,i)) = j;
    j=j+1;
  end
end

% identify boundary nodes in the global reordering and create gmapF and gmapB
%boundaryNodes = find(vmapM == vmapP);
nUnique = max(gmap);
boundaryNodes = unique(vmapB);
boundaryNodes = unique(gmap(boundaryNodes));
freeNodes = setdiff((1:nUnique)', boundaryNodes);
nFree = length(freeNodes);
nBoundary = length(boundaryNodes);
nTotal = nFree + nBoundary;

%[freeNodes; boundaryNodes] should be the final partitioning
perm = [freeNodes; boundaryNodes];
iperm(perm) = (1:nUnique);
gmap = iperm(gmap);
gmap = reshape(gmap, [Np,K])';

% create weighting vector to account for algebraic count of each node
%cmap = sum(gmap(:) == gmap(:)');
%cmap = reshape(cmap, [Np,K])';

% compute grid points in CG enumeration
xCG = zeros(nTotal,1); yCG = zeros(nTotal,1);
xCG(gmap') = x; yCG(gmap') = y;

end