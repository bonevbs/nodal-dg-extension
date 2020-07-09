function [bcx, bcy] = ElasticityIPDGbc2D(uxbc,uybc,qxbc,qybc)

% Purpose: Set up the discrete Poisson matrix directly
%          using LDG. The operator is set up in the weak form
Globals2D;

% build local face matrices
massEdge = zeros(Np,Np,Nfaces);
Fm = Fmask(:,1); faceR = r(Fm); V1D = Vandermonde1D(N, faceR);  massEdge(Fm,Fm,1) = inv(V1D*V1D');
Fm = Fmask(:,2); faceR = r(Fm); V1D = Vandermonde1D(N, faceR);  massEdge(Fm,Fm,2) = inv(V1D*V1D');
Fm = Fmask(:,3); faceS = s(Fm); V1D = Vandermonde1D(N, faceS);  massEdge(Fm,Fm,3) = inv(V1D*V1D');

% build DG right hand side
bcx = zeros(Np, K);
bcy = zeros(Np, K);

for k1=1:K
  if(~mod(k1,1000)) fprintf('ElasticityIPDGbc2D :: Processing element %d\n', k1), end;
  
  for f1=1:Nfaces
    
    if(BCType(k1,f1))
      
      Fm1 = Fmask(:,f1);
      fidM  = (k1-1)*Nfp*Nfaces + (f1-1)*Nfp + (1:Nfp)';
      
      id = 1+(f1-1)*Nfp + (k1-1)*Nfp*Nfaces;
      lnx = nx(id);  lny = ny(id); lsJ = sJ(id); hinv = Fscale(id);
      
      Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;
      Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;
      Dn1 = lnx*Dx  + lny*Dy ;
      
      mmE = lsJ*massEdge(:,:,f1);
      
      alpha = 100*2*(N+1)*(N+1)*hinv; % set penalty scaling
      switch(BCType(k1,f1))
        case {Dirichlet}
          bcx(:,k1) = bcx(:,k1) + alpha*mmE(:,Fm1)*uxbc(fidM);
          bcy(:,k1) = bcy(:,k1) + alpha*mmE(:,Fm1)*uybc(fidM);
        case {Neuman}
          bcx(:,k1) = bcx(:,k1) + 0.5*mmE(:,Fm1)*qxbc(fidM);
          bcy(:,k1) = bcy(:,k1) + 0.5*mmE(:,Fm1)*qybc(fidM);
        otherwise
      end
    end
  end
end
return
