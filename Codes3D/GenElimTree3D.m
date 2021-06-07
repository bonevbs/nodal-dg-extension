function elim_tree = GenElimTree3D(nmax, varargin)

% function elim_tree = GenElimTree3D(kmax, pmode)
% Purpose  : Given 3D mesh generate a nested dissection elimination tree;
%            This code is written from the bottom up and proceeds in three
%            steps:
%              1. Nested dissection is performed on the Elements.
%              2. Boundary elements are detected and tree is separated into
%                 boundaries and interior elements
%              3. Elements are replaced by corresponding nodes
% written by Boris Bonev

p = inputParser;
validPartitionModes = @(x) strcmp('alternating', x) || strcmp('vertical', x) || strcmp('horizontal', x) || strcmp('aspect', x);
validBoundaryTypes = @(x) strcmp('jagged', x) || strcmp('quads', x);
addRequired(p,'nmax', @(x) x>=0);
addParameter(p,'PartitionMode', 'alternating', validPartitionModes);
addParameter(p,'MaxLevels', 0, @(x) x>=0);
parse(p,nmax,varargin{:});

Globals3D;

% determine set of Elements and boundary elements
E = (1:K);
%EB = [];
EB = find(any(BCType,2)');

% generate hierarchy of elements
nd_tree = {};
nd_tree = nested_dissection_tree(nd_tree, E, EB, -1, 1, p.Results.nmax, p.Results.PartitionMode, p.Results.MaxLevels);
nd_tree = reverse_tree_indices(nd_tree);

% generate hierarchy of
elim_tree = gen_elimination_tree(nd_tree);
end

% build the tree recursively
function [nd_tree, i] = nested_dissection_tree(nd_tree, E, Bnd, i, level, nmax, pmode, maxlevels)
nd_tree = [{E; Bnd; i; [-1,-1]; level}, nd_tree];

Globals3D;

% handle children
if strcmp(pmode, 'horizontal')
  sep_mode='z';
elseif strcmp(pmode, 'vertical')
  sep_mode='x';
elseif strcmp(pmode, 'alternating')
  if rem(level,3) == 1
    sep_mode='x';
  elseif rem(level,3) == 2
    sep_mode='y';
  else
    sep_mode='z';
  end
elseif strcmp(pmode, 'aspect')
  xc = VX(EToV(E)); xc = xc(:); ax = abs(max(xc) - min(xc));
  yc = VY(EToV(E)); yc = yc(:); ay = abs(max(yc) - min(yc));
  zc = VZ(EToV(E)); zc = zc(:); az = abs(max(zc) - min(zc));
  if ax >= ay && ax >= az
    sep_mode='x';
  elseif ay >= ax && ay >= az
    sep_mode='y';
  else
    sep_mode='z';
  end
else
  error('unknown partition mode.')
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

Globals3D;

% if the critical number of Elements in the separator is reached just
% return all nodes as separator
if (size(EP, 2) < nmax)%|| pmax-pmin <= 1/16)
  EL = []; ER = []; BndL = []; BndR = [];
  return;
end

% decide on axis for separation
if strcmp(sep_mode, 'x')
  VP = VX;
elseif strcmp(sep_mode, 'y')
  VP = VY;
elseif strcmp(sep_mode, 'z')
  VP = VZ;
else
  error('unknown separation mode')
end

% find the bounding box and draw vertical line to separate
pmin = min(VP(EToV(EP))); pmax = max(VP(EToV(EP)));
psep = (pmax - pmin)/2 + pmin;

% find the elements which belong to the Separator and to Elements on
% the left/right
% VL1 = EP( VP(EToV(EP, 1)) < psep ); VR1 = EP( VP(EToV(EP, 1)) > psep );
% VL2 = EP( VP(EToV(EP, 2)) < psep ); VR2 = EP( VP(EToV(EP, 2)) > psep );
% VL3 = EP( VP(EToV(EP, 3)) < psep ); VR3 = EP( VP(EToV(EP, 3)) > psep );
% 
% EL = intersect( intersect( VL1, VL2, 'stable' ), VL3, 'stable' );
% ER = intersect( intersect( VR1, VR2, 'stable' ), VR3, 'stable' );
% ES = setdiff(setdiff(EP, EL, 'stable'), ER, 'stable');
% ES = union(ES, Bnd, 'stable');

%compute barycenter to figure out which elements go left/right
EL = EP( ( VP(EToV(EP, 1)) + VP(EToV(EP, 2)) + VP(EToV(EP, 3)) + VP(EToV(EP, 4)))/4 <= psep);
ER = EP( ( VP(EToV(EP, 1)) + VP(EToV(EP, 2)) + VP(EToV(EP, 3)) + VP(EToV(EP, 4)))/4 > psep);

% BndL = intersect(EL, ES, 'stable');
% BndR = intersect(ER, ES, 'stable');

BndL = EL(ismember(EToE(EL,1), ER) | ismember(EToE(EL,2), ER) | ismember(EToE(EL,3), ER) | ismember(EToE(EL,4), ER));
%BndR = intersect(unique(EToE(BndL,:), 'stable'), ER, 'stable')';
BndR = ER(ismember(EToE(ER,1), EL) | ismember(EToE(ER,2), EL) | ismember(EToE(ER,3), EL) | ismember(EToE(ER,4), EL));
BndL = union(BndL, intersect(Bnd, EL, 'stable'), 'stable');
BndR = union(BndR, intersect(Bnd, ER, 'stable'), 'stable');

% exit if one the sets ended up being just the boundary
if isempty(setdiff(EL,BndL)) || isempty(setdiff(ER,BndR))
  EL = []; ER = []; BndL = []; BndR = [];
end

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

Globals3D;

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
    EBnd = nd_tree{2,i};
    Bnd = ones([Np,1]) * (EBnd-1)*Np + (1:Np)' * ones([1,size(EBnd,2)]);
    Bnd = Bnd(:)';
    elim_tree{2,i} = Bnd;
    
    % construct the interior
    EInt = setdiff(nd_tree{1,i}, EBnd,'stable');
    Int = ones([Np,1]) * (EInt-1)*Np + (1:Np)' * ones([1,size(EInt,2)]);
    Int = setdiff(Int(:)',elim);
    elim = [elim,Int];
    elim_tree{1,i} = Int;
    
  end
end
end
