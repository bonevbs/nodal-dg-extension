% creates serialized datastructure for the elimination tree
function seq_tree = ElimTree2SequentialTree(elim_tree)

    nnodes = size(elim_tree,2);

    % extract tree structure
    fathers = [elim_tree{3,:}];
    sons    = [elim_tree{4,:}];
    lsons   = sons(1:2:2*nnodes);
    rsons   = sons(2:2:2*nnodes);
    
    % extract dimensions of index sets
    ninter  = cellfun(@length,{elim_tree{1,:}});
    nbound  = cellfun(@length,{elim_tree{2,:}});
    
    % extract index sets
    inter = zeros(max(ninter), nnodes);
    bound = zeros(max(nbound), nnodes);
    for i = 1:nnodes
        inter(1:ninter(i),i) = elim_tree{1,i}';
        bound(1:nbound(i),i) = elim_tree{2,i}';
    end

    seq_tree = struct('fathers', fathers, 'lsons', lsons, 'rsons', rsons, 'ninter', ninter, 'inter', inter, 'nbound', nbound, 'bound', bound);
end