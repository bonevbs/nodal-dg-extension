function elim_tree = GenElimTreeCG2D(nmax, varargin)

% function elim_tree = GenElimTree2DCG(kmax, pmode)
% Purpose  : Given 2D mesh generate a nested dissection elimination tree;
%            This code is written from the bottom up and proceeds in three
%            steps:
%              1. Nested dissection is performed on the Elements.
%              2. Boundary elements are detected and tree is separated into
%                 boundaries and interior elements
%              3. Elements are replaced by corresponding nodes
% written by Boris Bonev

p = inputParser;
validPartitionModes = @(x) strcmp('alternating', x) || strcmp('vertical', x) || strcmp('horizontal', x);
validBoundaryTypes = @(x) strcmp('jagged', x) || strcmp('quads', x);
addRequired(p,'nmax', @(x) x>=0);
addParameter(p,'PartitionMode', 'alternating', validPartitionModes);
addParameter(p,'MaxLevels', 0, @(x) x>=0);
parse(p,nmax,varargin{:});

Globals2D;
GlobalsCG2D;

E = (1:K);
%EB = [];
EB = any(BCType,2);
Bnd = unique(gmap(EB, :));
%Bnd = intersect(gmap(EB, :), 1:nFree, 'stable');

% generate hierarchy of elements
nd_tree = {};
nd_tree = nested_dissection_tree(nd_tree, E, Bnd, -1, 1, p.Results.nmax, p.Results.PartitionMode, p.Results.MaxLevels);
nd_tree = reverse_tree_indices(nd_tree);

% generate hierarchy of
elim_tree = gen_elimination_tree(nd_tree);
end

% build the tree recursively
function [nd_tree, i] = nested_dissection_tree(nd_tree, E, Bnd, i, level, nmax, pmode, maxlevels)
nd_tree = [{E; Bnd; i; [-1,-1]; level}, nd_tree];

% handle children
if strcmp(pmode, 'horizontal')
  sep_mode='y';
elseif strcmp(pmode, 'vertical')
  sep_mode='x';
else
  if rem(level,2); sep_mode='x'; else; sep_mode='y'; end
end
if level < maxlevels || maxlevels == 0
  [EL, ER, BndL, BndR] = separate(E, Bnd, sep_mode, nmax);
else
  EL = []; ER = []; BndL = []; BndR = [];
end
j = size(nd_tree,2)-1;

if ~isempty(ER)
  nd_tree{4, end-j}(2) = size(nd_tree,2);
  nd_tree = nested_dissection_tree(nd_tree, ER, BndR, j, level+1, nmax, pmode, maxlevels);
end
if ~isempty(EL)
  nd_tree{4, end-j}(1) = size(nd_tree,2);
  nd_tree = nested_dissection_tree(nd_tree, EL, BndL, j, level+1, nmax, pmode, maxlevels);
end
end

% function to find the separator nodes and two disjoint set o and returns 2 subtrees
function [EL, ER, BndL, BndR] = separate(EP, Bnd, sep_mode, nmax)

Globals2D;
GlobalsCG2D;

% if the critical number of Elements in the separator is reached just
% return all nodes as separator
if (size(EP, 2) < nmax)%|| pmax-pmin <= 1/16)
  EL = []; ER = []; BndL = []; BndR = [];
  return;
end

% decide on axis for separation
if strcmp(sep_mode, 'x')
  VP = VX;
  pos = xCG;
elseif strcmp(sep_mode, 'y')
  VP = VY;
  pos = yCG;
else
  error('unknown separation mode')
end

% find the bounding box and draw vertical line to separate
pmin = min(VP(EToV(EP))); pmax = max(VP(EToV(EP)));
psep = (pmax - pmin)/2 + pmin;

% find the elements which belong to the Separator and to Elements on
% the left/right
VL1 = EP( VP(EToV(EP, 1)) < psep ); VR1 = EP( VP(EToV(EP, 1)) > psep );
VL2 = EP( VP(EToV(EP, 2)) < psep ); VR2 = EP( VP(EToV(EP, 2)) > psep );
VL3 = EP( VP(EToV(EP, 3)) < psep ); VR3 = EP( VP(EToV(EP, 3)) > psep );

EL = intersect( intersect( VL1, VL2, 'stable' ), VL3, 'stable' );
ER = intersect( intersect( VR1, VR2, 'stable' ), VR3, 'stable' );
ES = setdiff(setdiff(EP, EL, 'stable'), ER, 'stable');

% get all DOFs on the separator and split them into two groups
Bnd = union(Bnd, gmap(ES,:),'stable');
indL = find(pos(Bnd) < psep);
indR = find(pos(Bnd) > psep);
indM = find(pos(Bnd) == psep);
indL = union(indL, indM(1:2:end));
indR = union(indR, indM(2:2:end));
BndL = Bnd(indL);
BndR = Bnd(indR);

% find the barycenters to have a conditioner for splitting 
%baryL = find( ( VP(EToV(EP, 1)) + VP(EToV(EP, 2)) + VP(EToV(EP, 3)) )/3 <= psep);
%baryR = find( ( VP(EToV(EP, 1)) + VP(EToV(EP, 2)) + VP(EToV(EP, 3)) )/3 > psep);


%BndR = intersect(Bnd, BndL, 'stable');

%EL = EP( ( VP(EToV(EP, 1)) + VP(EToV(EP, 2)) + VP(EToV(EP, 3)) )/3 <= psep);
%ER = EP( ( VP(EToV(EP, 1)) + VP(EToV(EP, 2)) + VP(EToV(EP, 3)) )/3 > psep);

end

% as indices where counted from the back of the tree fix that
function nd_tree = reverse_tree_indices(nd_tree)
nnodes = size(nd_tree, 2);
for i=1:nnodes
  if (nd_tree{3,i} >= 0); nd_tree{3,i} = nnodes-nd_tree{3,i}; end
  if (nd_tree{4,i}(1) >= 0); nd_tree{4,i}(1) = nnodes-nd_tree{4,i}(1); end
  if (nd_tree{4,i}(2) >= 0); nd_tree{4,i}(2) = nnodes-nd_tree{4,i}(2); end
end
end

% generation of the elimination tree
function elim_tree = gen_elimination_tree(nd_tree)

Globals2D;
GlobalsCG2D;

elim_tree = {};

levels = [nd_tree{5,:}];
L = max(levels);

% loop over all boxes from coarser to finer
elim = [];
for l = L:-1:1
  for i = 1:size(nd_tree,2)
    if (nd_tree{5,i} ~= l); continue; end
    
    %write current node into the tree
    elim_tree{3,i} = nd_tree{3,i};
    elim_tree{4,i} = nd_tree{4,i};
    elim_tree{5,i} = nd_tree{5,i};
    
    % Append boundary degrees of freedom
    Bnd = nd_tree{2,i};
    Bnd = Bnd(:)';
    Bnd = Bnd(Bnd <= nFree);
    elim_tree{2,i} = Bnd;
    
    % determine interior DOFs for elimination
    Int = setdiff(gmap(nd_tree{1,i},:), Bnd);
    Int = Int(:)';
    Int = Int(Int <= nFree);
    Int = setdiff(Int(:)',elim);
    elim = [elim,Int];
    elim_tree{1,i} = Int;
   
  end
end
end