function PlotSparsityReorder(A, sep_tree, tl)

    % generate permutation from separator tree
    [perm, ~] = ReorderSepTree(sep_tree, 'postorder');
    % get the dividers
    cellsz = cellfun(@(A)size(A,2),sep_tree,'uni',false);
    cellsz = cumsum([cellsz{1,:}]) + 1;
    levels = [sep_tree{4,:}];
    % select levels above d
    cellsz = cellsz(levels <= tl);
    %cellsz = unique(cellsz);
    
    % plot sparsity pattern and dividers
    figure();
    spy(A(perm,perm));
    grid on
    set(gca, 'GridLineStyle', '-')
    set(gca,'xtick',cellsz)
    set(gca,'ytick',cellsz)
    set(gca,'xticklabel',{[]}) 
    set(gca,'yticklabel',{[]})

end