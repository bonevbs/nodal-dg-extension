function PlotRanksReorder(A, seps, rkmode, rktol)

% function PlotRanksReorder(A, seps, rkmode, rktol)
% Purpose  : Given matrix A and seps, plot the ranks of all blocks induced by the partitioning;
%            seps is a cell array containing separators i.e. {[1,2,4], [3,5,6]}
% written by Boris Bonev 12/2018

    % reorder matrix according to separators
    perm = [seps{1,:}];
    % figure out how to scale the thing
    sepsz = cellfun(@(A)size(A,2),seps,'uni',false);
    sepsz = [sepsz{1,:}];
    nseps = size(sepsz, 2);

    % additional row/column required for pcolor plotting
    Aranks = zeros(nseps+1,nseps+1);

    for j=1:nseps
        for i=1:nseps
            rk = get_rank(A(seps{i}, seps{j}), rkmode, rktol);
            Aranks(i,j) = rk;
        end
    end

    %determine gridlines
    seploc = cumsum(sepsz);
    seploc = [0,seploc(1:end)];

    figure();
    pcolor(seploc,seploc,Aranks);
    colorbar()
    %colorbar(gray)
    axis ij
    axis square

    % add labels
    for jj = 1:nseps
        for ii = 1:nseps
            x = 0.5*(seploc(ii)+seploc(ii+1));
            y = 0.5*(seploc(jj)+seploc(jj+1));
            text(x, y, num2str(Aranks(ii,jj)), 'FontSize', 8, 'HorizontalAlignment', 'center');
        end
    end


end