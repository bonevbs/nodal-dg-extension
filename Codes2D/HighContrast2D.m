function [A,M] = HighContrast2D(a);

Globals2D;
g = zeros(K*Np,1);
A = spalloc(K*Np, K*Np, 3*Np);  M = spalloc(K*Np, K*Np, 3*Np); 

a = ones([g,Np,K]);

% Build matrix -- one column at a time
for i=1:K*Np
    g(i) = 1.0;
    gmat = reshape(g,Np,K);
    [Avec,Mvec] = HighContrastRHS2D(gmat, a);
   
    ids = find(Avec); A(ids,i) = Avec(ids);
    ids = find(Mvec); M(ids,i) = Mvec(ids);
    g(i)=0.0;
end
return

end