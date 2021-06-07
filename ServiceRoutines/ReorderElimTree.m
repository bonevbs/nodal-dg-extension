function reordered_tree = ReorderElimTree(elim_tree, perm)
    iperm(perm) = 1:length(perm);
    reordered_tree = {};
    assert(size(elim_tree,2) == size(perm,2),'Permutation vector and elimination tree need to have matching dimensions.');
    for i = 1:size(elim_tree,2)
        isep = perm(i);
        reordered_tree{1,i} = elim_tree{1,isep};
        reordered_tree{2,i} = elim_tree{2,isep};
        if (elim_tree{3,isep} ~= -1)
            reordered_tree{3,i} = iperm(elim_tree{3,isep});
        else
            reordered_tree{3,i} = -1;
        end
        reordered_tree{4,i} = [];
        for ison = elim_tree{4,isep}
            if (ison ~= -1)
                reordered_tree{4,i} = [reordered_tree{4,i}, iperm(ison)];
            else
                reordered_tree{4,i} = [reordered_tree{4,i}, -1];
            end
        end
        reordered_tree{5,i} = elim_tree{5,isep};
    end
end
