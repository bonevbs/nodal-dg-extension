function CheckElimTree(elim_tree, A)

% function CheckElimTree(elim_tree, A)
% Purpose  : Checks whether the elimination tree makes sense given the 
%            connectivity pattern of A.
% written by Boris Bonev 2019-09-17

assert(size(A,1) == size(A,2), 'A is not square.');

L = max([elim_tree{5,:}]);
N = size(A,2);
allDOF = 1:N;
remDOF = allDOF;
elmDOF = [];
nsep = size(elim_tree, 2);

correct = true;
for level=L:-1:1
    lvlDOF = [];
    for i=1:nsep
        if elim_tree{5,isep} ~= level
            continue
        else
            % check whether the interior only communicates with the 
            rest = setdiff(remDOF, [elim_tree{1,isep}, elim_tree{2,isep}]);
            
            lvlDOF = union(lvlDOF, elim_tree{1,isep});
        end
    end
end

end