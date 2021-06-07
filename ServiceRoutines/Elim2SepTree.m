function sep_tree = Elim2SepTree(elim_tree)

% function sep_tree = Elim2SepTree(elim_tree)
% Purpose  : converts elimination tree to separator tree.
%            this can be used with the old STRUMPACK code   
% written by Boris Bonev

sep_tree = elim_tree;

levels = [sep_tree{5,:}];
L = max(levels);
top_separator = union( [sep_tree{2,levels==1}], [sep_tree{1,levels==1}], 'stable');
sep_tree{1,levels==1} = top_separator;
% delete boundary DOFs
sep_tree(2,:) = [];

% just to make sure everything is alright
perm = [sep_tree{1,:}];
N = max(perm);
assert(all(ismember(1:N, perm)),'Some DOFs were lost in the conversion.')

end