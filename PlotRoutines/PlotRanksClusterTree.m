function PlotRanksClusterTree(A, cluster_tree, rkmode, rktol)

% function PlotRanksSepTree(A, cluster_tree, rkmode, rktol)
% Purpose  : Given matrix A and seps, plot the ranks of all blocks induced by the partitioning;
%            seps is a cell array containing separators i.e. {[1,2,4], [3,5,6]}
%            The clustering has to be contigious!
% written by Boris Bonev 12/2018

    maxL = max([cluster_tree{4,:}]);
    
    figure()
    % loop from top to bottom level and draw off-diagonal ranks
    for L=1:maxL
        % find leaves on current level
        leaves = find([cluster_tree{4,:}] == L);    
        for i=leaves
            % get children
            if (cluster_tree{3,i}(1) ~= -1); sepl = cluster_tree{1,cluster_tree{3,i}(1)}; else; sepl=[]; end
            if (cluster_tree{3,i}(2) ~= -1); sepr = cluster_tree{1,cluster_tree{3,i}(2)}; else; sepr=[]; end
            % first draw the off-diagonals if children exist
            if (~isempty(sepl) && ~isempty(sepr))
                % get children
                seploc = [sepl(1)-1, sepr(1)-1, sepr(end)];
                rk12 = get_rank(A(sepl, sepr), rkmode, rktol);
                rk21 = get_rank(A(sepr, sepl), rkmode, rktol);
                % do the plotting
                Aranks = [NaN, rk12, NaN; rk21, NaN, NaN; NaN, NaN, NaN];
                pcolor(seploc,seploc,Aranks);
                x = 0.5*(seploc(2)+seploc(3));
                y = 0.5*(seploc(1)+seploc(2));
                text(x, y, num2str(rk12), 'FontSize', 8, 'HorizontalAlignment', 'center');
                text(y, x, num2str(rk21), 'FontSize', 8, 'HorizontalAlignment', 'center');
                hold on
            % else determine the rank and draw diagonal block
            else
                if ~isempty(sepl)
                    sep = sepl;
                elseif ~isempty(sepr)
                    sep = sepr;
                else
                    sep = cluster_tree{1,i};
                end
                if isempty(sep); continue; end
                seploc = [sep(1)-1, sep(end)];
                rk = get_rank(A(sep, sep), rkmode, rktol);
                Aranks = [rk, NaN; NaN, NaN];
                pcolor(seploc,seploc,Aranks);
                x = 0.5*(seploc(1)+seploc(2));
                text(x, x, num2str(rk), 'FontSize', 8, 'HorizontalAlignment', 'center');
                hold on
            end
        end
    end
    
    colorbar()
    %colorbar(gray)
    axis ij
    axis square
    hold off
end