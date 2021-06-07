% reorders degrees of freedom in the elimination tree for systems of
% equations. Each degree of freedom is doubled, tripled etc depending on
% the dimension of the system.
function elim_tree = StackElimTree(elim_tree, N, dim)
  for isep = 1:size(elim_tree,2)
    int = elim_tree{1,isep};
    bnd = elim_tree{2,isep};
    int_stack = []; bnd_stack = [];
    for i=0:dim-1
      int_stack = [int_stack; (i*N)+int];
      bnd_stack = [bnd_stack; (i*N)+bnd];
    end
    elim_tree{1,isep} = int_stack(:)';
    elim_tree{2,isep} = bnd_stack(:)';
  end
end

% old version where DOFs are interleaved
% function elim_tree = StackElimTree(elim_tree, dim)
%   for i = 1:size(elim_tree,2)
%     int = elim_tree{1,isep};
%     bnd = elim_tree{2,isep};
%     int_stack = []; bnd_stack = [];
%     for i=(dim-1):-1:0
%       int_stack = [int_stack; dim*int-i];
%       bnd_stack = [bnd_stack; dim*int-i];
%     end
%     elim_tree{1,isep} = int_stack(:)';
%     elim_tree{2,isep} = bnd_stack(:)';
%   end
% end