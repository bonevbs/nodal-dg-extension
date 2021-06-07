function [perm, iperm] = ElimTree2Perm(elim_tree, mode)

    if strcmp(mode, 'postorder')
        perm = [elim_tree{1,:}, elim_tree{2,end}];
    elseif strcmp(mode, 'levels')
        % generate permutation for the elimination reordering
        levels = [elim_tree{5,:}];
        L = max(levels);
        perm = [];
        % go from top to bottom in the tree and arrange the reverse permutation
        for l=L:-1:1
            perm = [perm, elim_tree{1,levels==l}];
        end
        perm = [perm, elim_tree{2,levels==1}];
        % flip permutation
        %perm = fliplr(perm);
    else
        error('unknown reorder mode')
    end
    
    % generate the inverse permutation
    iperm(perm) = 1:length(perm);

end