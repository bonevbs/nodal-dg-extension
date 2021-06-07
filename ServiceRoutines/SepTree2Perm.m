function [perm, iperm] = ReorderSepTree(sep_tree, mode)


    if strcmp(mode, 'postorder')
        perm = [sep_tree{1,:}];
    elseif strcmp(mode, 'levels')
        % generate permutation for the elimination reordering
        levels = [sep_tree{4,:}];
        L = max(levels);
        perm = [];
        % go from top to bottom in the tree and arrange the reverse permutation
        for l=L:-1:1
            perm = [perm, sep_tree{1,levels==l}];
        end
        % flip permutation
        perm = fliplr(perm);
    else
        error('unknown reorder mode')
    end
    
    % generate the inverse permutation
    iperm(perm) = 1:length(perm);

end