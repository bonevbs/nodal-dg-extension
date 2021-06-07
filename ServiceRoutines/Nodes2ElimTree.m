function elim_tree = Nodes2ElimTree(NODES)
% function elim_tree = Nodes2ElimTree(NODES)
% Purpose  : converts NODES datastructure to elimination tree.
%            this can be used to convert Paolo's datastructures into mine.
% written by Boris Bonev

elim_tree = {};

% insert artificial top-level node
elim_tree{1,1} = [];
elim_tree{2,1} = [];
elim_tree{3,1} = -1;
elim_tree{4,1} = [];
elim_tree{5,1} = 1;

for isep=1:size(NODES,2)
    elim_tree{5,isep+1} = NODES{1,isep}+1;
    if (NODES{2,isep} == 0)
        elim_tree{4,1} = [elim_tree{4,1}, isep+1];
        elim_tree{3,isep+1} = 1;
        % append boundary nodes of sons to interior nodes of top-level
        % separator
        elim_tree{1,1} = union(elim_tree{1,1}, NODES{4,isep}, 'stable')';
    else
        elim_tree{3,isep+1} = NODES{2,isep}+1;
    end
    if (NODES{3,isep} == 0)
        elim_tree{4,isep+1} = [-1,-1];
    else
        elim_tree{4,isep+1} = NODES{3,isep}+1;
    end
    elim_tree{1,isep+1} = NODES{5,isep};
    elim_tree{2,isep+1} = NODES{4,isep};
end

perm = [];
perm = postorder_seps(perm, elim_tree, 1);
elim_tree = ReorderElimTree(elim_tree, perm);

end

% recursively permutation to order separators in postordering
function perm = postorder_seps(perm, elim_tree, isep)
    for ison = elim_tree{4,isep}
        if (ison ~= -1); perm = postorder_seps(perm, elim_tree, ison); end
    end
    perm = [perm, isep];
end