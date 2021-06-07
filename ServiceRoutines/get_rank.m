function rk = get_rank(B, rkmode, rktol)
    if strcmp(rkmode, 'rank')
        rk = rank(full(B), rktol);
    elseif strcmp(rkmode, 'relative rank')
        rk = rank(full(B), rktol)/min(size(B));
    elseif strcmp(rkmode, 'sprank')
        rk = sprank(B);
    elseif strcmp(rkmode, 'relative sprank')
        rk = sprank(B)/min(size(B));
    else
        error('unknown rank computation mode')
    end
end