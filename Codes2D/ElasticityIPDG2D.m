function [OPXX,OPXY,OPYX,OPYY,MM] = ElasticityIPDG2D(mu,lambda)

% Purpose: Set up the discrete Poisson matrix directly
%          using LDG. The operator is set up in the weak form
Globals2D;

% build local face matrices
massEdge = zeros(Np,Np,Nfaces);
Fm = Fmask(:,1); faceR = r(Fm);
V1D = Vandermonde1D(N, faceR);  massEdge(Fm,Fm,1) = inv(V1D*V1D');
Fm = Fmask(:,2); faceR = r(Fm);
V1D = Vandermonde1D(N, faceR);  massEdge(Fm,Fm,2) = inv(V1D*V1D');
Fm = Fmask(:,3); faceS = s(Fm);
V1D = Vandermonde1D(N, faceS);  massEdge(Fm,Fm,3) = inv(V1D*V1D');

% build local volume mass matrix
MassMatrix = invV'*invV;

% build DG derivative matrices
MM  = zeros(K*Np*Np, 3);
OPXX = zeros(K*Np*Np*(1+Nfaces), 3); OPXY = zeros(K*Np*Np*(1+Nfaces), 3);
OPYX = zeros(K*Np*Np*(1+Nfaces), 3); OPYY = zeros(K*Np*Np*(1+Nfaces), 3);

% global node numbering
entries = (1:Np*Np)'; entriesMM = (1:Np*Np)';
for k1=1:K
  if(~mod(k1,1000)) fprintf('ElasticityIPDG2D :: Processing element %d\n', k1), end;
  rows1 = ((k1-1)*Np+1:k1*Np)'*ones(1,Np); cols1 = rows1';
  
  % Build local operators
  Dx = rx(1,k1)*Dr + sx(1,k1)*Ds; Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;
  
  % elasticity parameters
  Mu1 = diag(mu(:,k1)); Lambda1 = diag(lambda(:,k1));
  
  % interior variational terms (Term 1)
  ODXX = J(1,k1)*(Dx'*MassMatrix*(2*Mu1 + Lambda1)*Dx + Dy'*MassMatrix*Mu1*Dy);
  ODYY = J(1,k1)*(Dy'*MassMatrix*(2*Mu1 + Lambda1)*Dy + Dx'*MassMatrix*Mu1*Dx);
  ODXY = J(1,k1)*(Dx'*MassMatrix*Lambda1*Dy + Dy'*MassMatrix*Mu1*Dx);
  ODYX = J(1,k1)*(Dy'*MassMatrix*Lambda1*Dx + Dx'*MassMatrix*Mu1*Dy);
  
  % Build element-to-element parts of operator
  for f1=1:Nfaces
    k2 = EToE(k1,f1); f2 = EToF(k1,f1);
    % elasticity parameters
    Mu2 = diag(mu(:,k2)); Lambda2 = diag(lambda(:,k2));
    
    rows2 = ((k2-1)*Np+1:k2*Np)'*ones(1,Np); cols2 = rows2';
    
    fidM  = (k1-1)*Nfp*Nfaces + (f1-1)*Nfp + (1:Nfp);
    vidM = vmapM(fidM); Fm1 = mod(vidM-1,Np)+1;
    vidP = vmapP(fidM); Fm2 = mod(vidP-1,Np)+1;
    
    id = 1+(f1-1)*Nfp + (k1-1)*Nfp*Nfaces;
    lnx = nx(id); lny = ny(id); lsJ = sJ(id);
    hinv = max(Fscale(id), Fscale(1+(f2-1)*Nfp, k2));
    
    Dx2 = rx(1,k2)*Dr + sx(1,k2)*Ds; Dy2 = ry(1,k2)*Dr + sy(1,k2)*Ds;
    
    %Dn1 = lnx*Dx  + lny*Dy ;
    %Dn2 = lnx*Dx2 + lny*Dy2;
    
    Dnxx1 = lnx*Dx; Dnxy1 = lnx*Dy; Dnyx1 = lny*Dx; Dnyy1 = lny*Dy;
    Dnxx2 = lnx*Dx2; Dnxy2 = lnx*Dy2; Dnyx2 = lny*Dx2; Dnyy2 = lny*Dy2;
    
    mmE = lsJ*massEdge(:,:,f1);
    
    % linear elasticity
    % interior variational terms
    
    % matrix for off-diagonal entries
    OOXX = zeros(Np); OOXY = zeros(Np); OOYX = zeros(Np); OOYY = zeros(Np);
    
    alpha = 10000*2*(N+1)*(N+1)*hinv;
    %alphaN = 100*2*(N+1)*(N+1)*hinv;
    %alphaT = 100*2*(N+1)*(N+1)*hinv;
    switch(BCType(k1,f1))
      case {Dirichlet}
        
        % sigmaN is not defined on Dirichlet boundaries
        ODXX = ODXX - mmE*(2*Mu1*Dnxx1 + Mu1*Dnyy1 + Lambda1*Dnxx1);
        ODXY = ODXY - mmE*(Mu1*Dnyx1 + Lambda1*Dnxy1);
        ODYX = ODYX - mmE*(Mu1*Dnxy1 + Lambda1*Dnyx1);
        ODYY = ODYY - mmE*(Mu1*Dnxx1 + 2*Mu1*Dnyy1 + Lambda1*Dnyy1);

        % penalization Term 4
        ODXX = ODXX + alpha*mmE;
        %ODXY = ODXY + alpha*mmE;
        %ODYX = ODYX + alpha*mmE;
        ODYY = ODYY + alpha*mmE;

      case {Neuman} 

      otherwise
        
        % diagonal contributions of Term 2
        ODXX = ODXX - 0.5*mmE*(2*Mu1*Dnxx1 + Mu1*Dnyy1 + Lambda1*Dnxx1);
        ODXY = ODXY - 0.5*mmE*(Mu1*Dnyx1 + Lambda1*Dnxy1);
        ODYX = ODYX - 0.5*mmE*(Mu1*Dnxy1 + Lambda1*Dnyx1);
        ODYY = ODYY - 0.5*mmE*(Mu1*Dnxx1 + 2*Mu1*Dnyy1 + Lambda1*Dnyy1);
        
        % off-diagonal contribution of Term 2
        OOXX(Fm1,:) = OOXX(Fm1,:) - 0.5*mmE(Fm1,Fm1)*(2*Mu2(Fm2,Fm2)*Dnxx2(Fm2,:) + Mu2(Fm2,Fm2)*Dnyy2(Fm2,:) + Lambda2(Fm2,Fm2)*Dnxx2(Fm2,:));
        OOXY(Fm1,:) = OOXY(Fm1,:) - 0.5*mmE(Fm1,Fm1)*(Mu2(Fm2,Fm2)*Dnyx2(Fm2,:) + Lambda2(Fm2,Fm2)*Dnxy2(Fm2,:));
        OOYX(Fm1,:) = OOYX(Fm1,:) - 0.5*mmE(Fm1,Fm1)*(Mu2(Fm2,Fm2)*Dnxy2(Fm2,:) + Lambda2(Fm2,Fm2)*Dnyx2(Fm2,:));
        OOYY(Fm1,:) = OOYY(Fm1,:) - 0.5*mmE(Fm1,Fm1)*(Mu2(Fm2,Fm2)*Dnxx2(Fm2,:) + 2*Mu2(Fm2,Fm2)*Dnyy2(Fm2,:) + Lambda2(Fm2,Fm2)*Dnyy2(Fm2,:));
        
        % penalization Term 4
        ODXX = ODXX + alpha*mmE;
        %ODXY = ODXY + alpha*mmE;
        %ODYX = ODYX + alpha*mmE;
        ODYY = ODYY + alpha*mmE;
        
        % off-diagonal penalization Term 4
        OOXX(:,Fm2) = OOXX(:,Fm2) - alpha*mmE(:,Fm1);
        %OOXY(:,Fm2) = OOXY(:,Fm2) - alpha*mmE(:,Fm1);
        %OOYX(:,Fm2) = OOYX(:,Fm2) - alpha*mmE(:,Fm1);
        OOYY(:,Fm2) = OOYY(:,Fm2) - alpha*mmE(:,Fm1);
        
%         % diagonal contributions of normal part of Term 4
%         ODXX = ODXX + alphaN*lnx*mmE*lnx;
%         ODXY = ODXY + alphaN*lnx*mmE*lnx;
%         ODYX = ODYX + alphaN*lnx*mmE*lnx;
%         ODYY = ODYY + alphaN*lnx*mmE*lnx;
% 
%         % off-diagonal contribution of normal part of Term 4 (only at internal and Dirichlet edges)
%         OOXX(:,Fm2) = OOXX(:,Fm2) - alphaN*lnx*mmE(:,Fm1)*lnx;
%         OOXY(:,Fm2) = OOXY(:,Fm2) - alphaN*lnx*mmE(:,Fm1)*lnx;
%         OOYX(:,Fm2) = OOYX(:,Fm2) - alphaN*lnx*mmE(:,Fm1)*lnx;
%         OOYY(:,Fm2) = OOYY(:,Fm2) - alphaN*lnx*mmE(:,Fm1)*lnx;
        
        % move entries to the appropriate place
        OPXX(entries(:), :) = [rows1(:), cols2(:), OOXX(:)];
        OPXY(entries(:), :) = [rows1(:), cols2(:), OOXY(:)];
        OPYX(entries(:), :) = [rows1(:), cols2(:), OOYX(:)];
        OPYY(entries(:), :) = [rows1(:), cols2(:), OOYY(:)];
        entries = entries + Np*Np;
    end
  end
  OPXX(entries(:), :) = [rows1(:), cols1(:), ODXX(:)];
  OPXY(entries(:), :) = [rows1(:), cols1(:), ODXY(:)];
  OPYX(entries(:), :) = [rows1(:), cols1(:), ODYX(:)];
  OPYY(entries(:), :) = [rows1(:), cols1(:), ODYY(:)];
  MM(entriesMM(:), :) = [rows1(:), cols1(:), J(1,k1)*MassMatrix(:)];
  entries = entries + Np*Np; entriesMM = entriesMM + Np*Np;
end

OPXX = OPXX(1:max(entries)  -Np*Np,:);  OPXX = myspconvert(OPXX, Np*K, Np*K, 1e-15);
OPXY = OPXY(1:max(entries)  -Np*Np,:);  OPXY = myspconvert(OPXY, Np*K, Np*K, 1e-15);
OPYX = OPYX(1:max(entries)  -Np*Np,:);  OPYX = myspconvert(OPYX, Np*K, Np*K, 1e-15);
OPYY = OPYY(1:max(entries)  -Np*Np,:);  OPYY = myspconvert(OPYY, Np*K, Np*K, 1e-15);
MM   = MM(1:max(entriesMM)-Np*Np,:);    MM   = myspconvert(MM, Np*K, Np*K, 1e-15);
return
