function sep_tree = GenNestedDissection3D(nmax)

% function sep_tree = GenNestedDissection2D()
% Purpose  : Given 3D mesh generate a nested dissection separator tree;   
% written by Boris Bonev

    Globals2D;

    sep_tree = {};
    sep_tree = gen_sep_tree(sep_tree, (1:K)', -1, 1, nmax);
    sep_tree = reverse_indices(sep_tree);

end

% build the tree recursively
function [SepTree, i] = gen_sep_tree(SepTree, Sep, i, level, nmax)
    if rem(level,3) == 0; sep_mode='z'; elseif rem(level,3) == 1; sep_mode='y'; else; sep_mode='x'; end
    [Sep, EL, ER] = separate(Sep, sep_mode, nmax);
    % append cleaned up separator
    %remove if
    if size(SepTree,2) > 0; Sep = setdiff(Sep, [SepTree{1, :}], 'stable'); end
    SepTree = [{Sep'; i; [-1,-1]; level}, SepTree];
    j = size(SepTree,2)-1;
    
    if ~isempty(ER)
        SepTree{3, end-j}(2) = size(SepTree,2);
        SepTree = gen_sep_tree(SepTree, ER, j, level+1, nmax);
    end
    if ~isempty(EL)
        SepTree{3, end-j}(1) = size(SepTree,2);
        SepTree = gen_sep_tree(SepTree, EL, j, level+1, nmax);
    end
end

% function to find the separator nodes and two disjoint set o and returns 2 subtrees
function [Sep, EL, ER] = separate(EP, sep_mode, nmax)

    Globals3D;
    
    % decide on axis for separation
    if strcmp(sep_mode, 'x')
        VP = VX;
    elseif strcmp(sep_mode, 'y')
        VP = VY;
    elseif strcmp(sep_mode, 'z')
        VP = VZ;
    else
        error('unknown separation mode')
    end
    
    % find the bounding box and draw vertical line to separate
    pmin = min(VP(EToV(EP))); pmax = max(VP(EToV(EP)));
    psep = (pmax - pmin)/2 + pmin;
    
    %compute barycenter to figure out which elements go left/right
    EL = find( ( VP(EToV(EP, 1)) + VP(EToV(EP, 2)) + VP(EToV(EP, 3)) + VP(EToV(EP, 4)))/4 <= psep)';
    ER = find( ( VP(EToV(EP, 1)) + VP(EToV(EP, 2)) + VP(EToV(EP, 3)) + VP(EToV(EP, 4)))/4 > psep)';
    EL = EP(EL); ER = EP(ER);

    % if the critical number of Elements in the separator is reached just
    % return all nodes as separator
    if (size(EL, 1)*Nfp < nmax || size(ER, 1)*Nfp < nmax )%|| pmax-pmin <= 1/16)
        EL = []; ER = [];
        Sep = ones([Np,1]) * (EP-1)'*Np + (1:Np)' * ones([1,size(EP,1)]);
        Sep = Sep(:);
        return
        
    % if separation is to be performed, draw a separator and split elements
    % into two disjoint sets. find the overlapping nodes which lie on the
    % interfaces
    else
        
        % check El and Er are disjoint
        assert(isempty(intersect(EL,ER)),'El and Er are not disjoint.')
        
        % select indices corresponding to nodes on the interfaces of El/Er
        NdL = ones([Nfp*Nfaces,1]) * (EL-1)'*Nfp*Nfaces + (1:(Nfp*Nfaces))' * ones([1,size(EL,1)]);
        NdR = ones([Nfp*Nfaces,1]) * (ER-1)'*Nfp*Nfaces + (1:(Nfp*Nfaces))' * ones([1,size(ER,1)]);
        NdL = NdL(:);
        NdR = NdR(:);

        Sep = union( intersect(vmapP(NdL),vmapM(NdR),'stable'), intersect(vmapM(NdL),vmapP(NdR),'stable') , 'stable');
        
        return
    end
        
end

% as indices where counted from the back of the tree traverse and fix
% indices
function SepTree = reverse_indices(SepTree)
    nbsep = size(SepTree, 2);
    for i=1:nbsep
        if (SepTree{2,i} >= 0); SepTree{2,i} = nbsep-SepTree{2,i}; end
        if (SepTree{3,i}(1) >= 0); SepTree{3,i}(1) = nbsep-SepTree{3,i}(1); end
        if (SepTree{3,i}(2) >= 0); SepTree{3,i}(2) = nbsep-SepTree{3,i}(2); end
    end
end
