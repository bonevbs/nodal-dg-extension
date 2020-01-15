function [A,M] = PoissonCG2D()

% function [A,M] = PoissonCG2D()
% Purpose: Set up matrix for 2D Poisson equation based on stabilized 
%          internal fluxes on symmetric form

Globals2D;
GlobalsCG;

A = sparse(nTotal, nTotal);  M = sparse(nTotal, nTotal); 

% Build matrices -- local matrices in DG formulation are written into the
% corresponding place in CG enumeration
for k=1:K
  Sx = (diag(rx(:,k))*Dr)'*MassMatrix*(diag(rx(:,k))*Dr)...
     + (diag(rx(:,k))*Dr)'*MassMatrix*(diag(sx(:,k))*Ds)...
     + (diag(sx(:,k))*Ds)'*MassMatrix*(diag(rx(:,k))*Dr)...
     + (diag(sx(:,k))*Ds)'*MassMatrix*(diag(sx(:,k))*Ds);
  Sy = (diag(ry(:,k))*Dr)'*MassMatrix*(diag(ry(:,k))*Dr)...
     + (diag(ry(:,k))*Dr)'*MassMatrix*(diag(sy(:,k))*Ds)...
     + (diag(sy(:,k))*Ds)'*MassMatrix*(diag(ry(:,k))*Dr)...
     + (diag(sy(:,k))*Ds)'*MassMatrix*(diag(sy(:,k))*Ds);
  
  %J(1,k)*(Sx+Sy)
  M(gmap(k,:), gmap(k,:)) = M(gmap(k,:), gmap(k,:)) + J(1,k)*MassMatrix;
  A(gmap(k,:), gmap(k,:)) = A(gmap(k,:), gmap(k,:)) + J(1,k)*(Sx+Sy);
end
return