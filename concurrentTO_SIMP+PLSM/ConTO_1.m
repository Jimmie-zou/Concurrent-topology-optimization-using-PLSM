function ConTO_plsm %(nx,ny,vol_con,er,rmin)
scale = 1;
nx = 120*scale; ny = 50*scale;
nelx = 40; nely = 40;
EW=1; EH=1; cIMQ=1/(EW*EH/nelx/nely);
DomainSize = [EW EH];
EleNumPerRow=nelx;  EleNumPerCol=nely;
E0 = 2.01e7; nu = 0.3; 
er = 0.04; rmin = 5*scale^1.5; 
vol_ma = 0.6; vol_mi = 0.5;
%% INITIALIZE PATTERN
x = ones(ny,nx); vol = 1; iter = 0; change = 1; c = [];
% initial PUC and \alpha
HoleRadius = 9;
rmin=2.5;
passive = zeros(nely+1,nelx+1);
passive(1:2,1:2) = 0;        passive(1:2,end-1:end) = 0;
passive(end-1:end,1:2) = 0;  passive(end-1:end,end-1:end) = 0;
aspect = DomainSize(2)/EleNumPerCol;

[Phi0, RBFGridXYi, FEANodeXYi] = IniDesign_New2(nelx, nely, DomainSize, HoleRadius);
NodeXYi = FEANodeXYi;  % all element nodes' coordinates: n-by-2 dimension
CentXYi = RBFGridXYi;    % all RBF knots' coordinates: n-by-2 dimension
[A] = CollecMatrix_IMQ(NodeXYi, CentXYi, cIMQ, DomainSize, nelx, nely);
% save PLSM A
% load('PLSM.mat', 'A')
Phi0 = reshape(Phi0,nely+1,nelx+1); Phi0(find(passive)) = 1/3*max(max(Phi0));
Phi0 = reshape(Phi0,(nely+1)*(nelx+1),1);
Phi = Phi0; 
f1 = Phi0;                 % the initial value "f" over knotes : (nely+1)*(nelx+1)-by-1
[R1,R2] = dwt(f1,'db1');
fx=[R1;R2]; Bx = A\fx; 
coef = idwt(Bx(1:length(R1)),Bx((length(R1)+1):end),'db1');
coeff1 = coef(1:(nelx+1)*(nely+1));
Coef= coeff1;
%% INDEXING NODES AND ELEMENT
nodenrs = reshape(1:(1+nx)*(1+ny),1+ny,1+nx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nx*ny,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*ny+[2 3 0 1] -2 -1],nx*ny,1);
% % stiffness matrix of solid material
% lambda = nu*E0/(1+nu)/(1-2*nu);
% mu = E0/2/(1+nu);    lambda = 2*lambda*mu/(lambda+2*mu);
% C0 = lambda*[1 1 0; 1 1 0; 0 0 0]+mu*[2 0 0; 0 2 0; 0 0 1];
% [KE] = StiffnessMatrix(C0,EW,EH);
iK = reshape(kron(edofMat,ones(8,1))',64*nx*ny,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nx*ny,1);
%% DEFINE LOADS AND SUPPORTS (CANTILEVER)
F = sparse(2*(nx+1)*(ny+1),1,-1e4,2*(ny+1)*(nx+1),1);
fixeddofs = 1:2*(ny+1);
% F = sparse(2*(ny+1)*(nx/2)+2,1,-1,2*(ny+1)*(nx+1),1);
% fixeddofs = [2*(ny+1),2*(ny+1)-1,2*(ny+1)*(nx+1)];
% F = sparse(2,1,-1e4,2*(ny+1)*(nx+1),1);
% fixeddofs = [1:2:2*(ny+1),2*(ny+1)*(nx+1)];
U = zeros(2*(ny+1)*(nx+1),1);
alldofs = 1:2*(ny+1)*(nx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
[H,Hs]=compute_filter(nx,ny,7.5);
%%
color=[255 255 255;     0 0 239;     0 0 255;     0 16 255;    0 32 255;   0 48 255;   0 64 255;   0 80 255;    0 96 255;  
        0 111 255;   0 128 255;   0 143 255;   0 159 255;   0 175 255;  0 191 255;  0 207 255;  0 223 255;   0 239 255;   
        0 255 255;   16 255 239;  32 255 223;  48 255 207;  64 255 191; 80 255 175; 96 255 159; 111 255 143; 128 255 128; 
        143 255 111; 159 255 96;  175 255 80;  191 255 64;  207 255 48; 223 255 32; 239 255 16; 255 239 0;   255 223 0;   
        255 207 0;   255 191 0;   255 175 0;   255 159 0;   255 143 0;  255 128 0;  255 111 0;  255 96 0;    255 80 0;   
        255 64 0;    255 48 0;    255 32 0;    255 16 0;    255 0 0;    239 0 0;    223 0 0;    207 0 0;     191 0 0;]./255;
%% make dir
str1 = strcat('Micro_Results_',num2str(EleNumPerRow),'x',num2str(EleNumPerCol),'_', num2str(vol_mi));
parent_dir_name =strcat('Experiment_results_', num2str(vol_ma),'\');
Micro_result_dir = strcat(parent_dir_name,str1);
mkdir(Micro_result_dir);
str2 = strcat('Macrostructure',num2str(nx),'x',num2str(ny),'_v',num2str(vol_ma));
Macro_results_dir = strcat(parent_dir_name,str2);
mkdir(Macro_results_dir);
%% MAIN OPTIMIZATION LOOP
maxit = 150;
while change > 0.001 && iter < maxit
    iter = iter + 1; 
    [dCH11,dCH12,dCH13,dCH22,dCH23,dCH33,Ke,A0,CH] = FEACalculation_CM(Phi,nelx,nely,EW,EH,EW,EH,E0,nu);
    vol = max(vol*(1-er),vol_ma);
    if iter >1; olddc = dc; end
    % FINITE ELEMENT ANALYSIS
    xmin = 1e-6;
    x_tilde = max(x, xmin);
    sK = reshape(Ke(:)*x_tilde(:)',64*nx*ny,1); 
    K = sparse(iK,jK,sK); K = (K+K')/2;
    % U(freedofs) = decomposition(K(freedofs,freedofs),'chol','lower')\F(freedofs);  
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    ce = reshape(sum((U(edofMat)*Ke).*U(edofMat),2),ny,nx);
    c(iter) = sum(sum(x.*ce))/2;
%     c = [c, 0.5.*sum(sum(x_sandwich*E0.*ce))];
    % MACROSCALE SENSITIVITY ANALYSIS
    dc = -x.*ce;                                                                                                                 
    dv = ones(ny,nx);
    dc(:) = H*dc(:)./Hs;
    if iter > 1; dc = (dc+olddc)/2.; end  % STABILIZATION OF EVOLUTIONARY PROCESS 
    % updating
    x=ADDDLE(ny,nx,dc,dv,vol);

    [~,Phi,Coef] = TopMMALSM(EleNumPerRow, EleNumPerCol,CentXYi,iter, DomainSize,vol_mi, passive,...
        nelx, nely, Coef, Phi, A, rmin, aspect,cIMQ,A0,reshape(Ke(:),8,8),U,nx,ny,dCH11,dCH12,dCH13,dCH22,dCH23,dCH33,EW,EH,CH,edofMat,x,Micro_result_dir);

    % PRINT RESULTS AND PLOT DENSITIES
    if iter > 15
        % change = abs(sum(c(iter-9:iter-5))-sum(c(iter-4:iter)))/sum(c(iter-4:iter));
        change = abs(sum(c(iter-7:iter-4))-sum(c(iter-3:iter)))/sum(c(iter-3:iter));
    end
    disp([' It.: ' sprintf('%2i',iter) ' Obj.: ' sprintf('%6.4f',c(iter)) ' Vma.: ' sprintf('%6.4f',sum(x(:))/nx/ny) ' Vmi.:' sprintf('%6.4f',A0) ' ch.: ' sprintf('%6.4f',change)])
    h1=figure(101);clf; colormap(gray); imagesc(1-x); clim([0 1]);axis equal tight off; pause(1e-6);
    FileName=[Macro_results_dir,'\Macro',int2str(iter),'.tif'];
    saveas(h1,FileName);
%     FileName=[parent_dir_name,'\Obj.txt'];
%     save(FileName, 'Obj','-ascii');
%     VVol(itte,:) = fval+volfrac;
%     FileName=[parent_dir_name,'\VVol.txt'];
%     save(FileName, 'VVol','-ascii');
end
end
%% PREPARE FILTER
function [H,Hs]=compute_filter(nelx,nely,rmin)
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));  sH = zeros(size(iH));
index=0.;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        [i2,j2] = ndgrid(max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx),max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely));
        e2 = (i2(:)-1)*nely+j2(:);
        iH(index + (1:numel(e2))) = e1;
        jH(index + (1:numel(e2))) = e2;
        sH(index + (1:numel(e2))) = max(0,rmin-sqrt((i1-i2(:)).^2+(j1-j2(:)).^2));
        index = index + numel(e2);
    end
end
H = sparse(iH,jH,sH); Hs = sum(H,2);
end
%% UPDATING
function [x] = ADDDLE(nely,nelx,dc,dv,vol)
    l1 = min(min(-dc./dv)); l2 = max(max(-dc./dv));
    x = max(0,sign(-dc./dv-l1));
    while ((l2-l1)/l2 > 1e-9)
        th = (l1+l2)/2;
        x = max(0,sign(-dc./dv-th));
        if sum(sum(x))-vol*(nelx*nely) > 0
            l1 = th;
        else
            l2 = th;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%