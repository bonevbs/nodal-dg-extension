function [cluster_tree, perm, iperm] = Sep2ClusterTree(sep_tree)

    [perm, iperm] = ReorderSepTree(sep_tree, 'postorder');
    
    %iperm([sep_tree{1,:}])
    cluster_tree = {};
    cl_start = 1;
    cl_end   = size(sep_tree,2);
    cluster_tree = gen_cluster_tree(cluster_tree, sep_tree, iperm, cl_start, cl_end, -1, 1);
    
end


function cluster_tree = gen_cluster_tree(cluster_tree, sep_tree, iperm, cl_start, cl_end, iparent, level)

% this function takes the current separator and appends the respective
% cluster tree in the form:
%
%             |
%            cl
%           /  \
%          /    \
%      sep_comp sep
%      /     \
%     /       \
%   cll       clr

    % take care of current cluster
    i = size(cluster_tree,2)+1;
    cl = [sep_tree{1,cl_start:cl_end}];
    cl = iperm(cl);
    cluster_tree{1,i} = cl;
    cluster_tree{2,i} = iparent;
    cluster_tree{4,i} = level;
    cluster_tree{3,i} = [-1, -1];
    if isempty([sep_tree{1,cl_start:cl_end-1}]); return; end
    cluster_tree{3,i} = [i+1, i+2];
    
    sep_comp = [sep_tree{1,cl_start:cl_end-1}];
    sep_comp = iperm(sep_comp);
    cluster_tree{1,i+1} = sep_comp;
    cluster_tree{2,i+1} = i;
    cluster_tree{3,i+1} = [-1,-1];
    cluster_tree{4,i+1} = level+1;
    
    sep = sep_tree{1,cl_end};
    sep = iperm(sep);
    cluster_tree{1,i+2} = sep;
    cluster_tree{2,i+2} = i;
    cluster_tree{3,i+2} = [-1,-1];
    cluster_tree{4,i+2} = level+1;
    
    % get left and right child of sep to find the appropriate range
    chl = sep_tree{3,cl_end}(1);
    chr = sep_tree{3,cl_end}(2);
    cll_end = 0;
    if chl ~= -1
        cll_start = cl_start;
        cll_end   = chl;
        cluster_tree{3,i+1}(1) = size(cluster_tree,2)+1;
        cluster_tree = gen_cluster_tree(cluster_tree, sep_tree, iperm, cll_start, cll_end, i+1, level+2);
    end
    if chr ~= -1
        clr_start = chl+1;
        clr_end   = chr;
        cluster_tree{3,i+1}(2) = size(cluster_tree,2)+1;
        cluster_tree = gen_cluster_tree(cluster_tree, sep_tree, iperm, clr_start, clr_end, i+1, level+2);
    end
end