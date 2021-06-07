function cluster_tree = GenClusterTreeBisection(seps)

% function cluster_tree = GenClusterTreeBisection(seps)
% Purpose  : Given a list of separators generate a cluster tree using bisection;   
% written by Boris Bonev 12/2018

% same binary tree datastructure as in the rest of the code
% cluster_tree{1,i}: contains ith index set
% cluster_tree{2,i}: contains a pointer to the parent of i-th node
% cluster_tree{3,i}: contains pointers to the children of i-th node
% cluster_tree{4,i}: contains the level of the current node

    % initialize cluster tree
    cluster_tree = {};
    % build cluster tree recursively
    cluster_tree = gen_cluster_tree(cluster_tree, seps, -1, 1);
    cluster_tree = reverse_indices(cluster_tree);
end

function cluster_tree = gen_cluster_tree(cluster_tree, seps, i, level)
    cluster_tree = [{[seps{1,:}]; i; [-1,-1]; level}, cluster_tree];
    j = size(cluster_tree,2)-1;   
    mid = ceil(size(seps,2)/2);
    if mid == size(seps,2)
        return
    end
    sl = {seps{1:mid}};
    sr = {seps{mid+1:end}};
    if ~isempty(sr)
        cluster_tree{3, end-j}(2) = size(cluster_tree,2);
        cluster_tree = gen_cluster_tree(cluster_tree, sr, j, level+1);
    end
    if ~isempty(sl)
        cluster_tree{3, end-j}(1) = size(cluster_tree,2);
        cluster_tree = gen_cluster_tree(cluster_tree, sl, j, level+1);
    end
end


% as indices where counted from the back of the tree traverse and fix
% indices
function cluster_tree = reverse_indices(cluster_tree)
    nbsep = size(cluster_tree, 2);
    for i=1:nbsep
        if (cluster_tree{2,i} >= 0); cluster_tree{2,i} = nbsep-cluster_tree{2,i}; end
        if (cluster_tree{3,i}(1) >= 0); cluster_tree{3,i}(1) = nbsep-cluster_tree{3,i}(1); end
        if (cluster_tree{3,i}(2) >= 0); cluster_tree{3,i}(2) = nbsep-cluster_tree{3,i}(2); end
    end
end