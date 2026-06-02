%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function top88(nx,ny,volfrac,penal,rmin,ft)
%% MATERIAL PROPERTIES
nx = 60; ny=30; volfrac = 0.5; 
nelx = 40; nely = 40;
penal = 3; rmin= 2.5; ft = 1;
EW=1; EH=1; cIMQ=1/(EW*EH/nelx/nely);
DomainSize = [EW EH];
EleNumPerRow=nelx;  EleNumPerCol=nely;
E0 = 2.01e7; nu = 0.3; 
vol_ma = 0.6; vol_mi = 0.5;
%% PREPARE FINITE ELEMENT ANALYSIS
nodenrs = reshape(1:(1+nx)*(1+ny),1+ny,1+nx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nx*ny,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*ny+[2 3 0 1] -2 -1],nx*ny,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nx*ny,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nx*ny,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
% F = sparse(2,1,-1e4,2*(ny+1)*(nx+1),1);
% fixeddofs = union([1:2:2*(ny+1)],[2*(nx+1)*(ny+1)]);
F = sparse(2*(ny+1)*(nx)+ny+2,1,-1e4,2*(ny+1)*(nx+1),1);
fixeddofs = 1:2*(ny+1);
U = zeros(2*(ny+1)*(nx+1),1);
alldofs = [1:2*(ny+1)*(nx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nx*ny*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));  sH = zeros(size(iH));
index=0.;
for i1 = 1:nx
    for j1 = 1:ny
        e1 = (i1-1)*ny+j1;
        [i2,j2] = ndgrid(max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nx),max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),ny));
        e2 = (i2(:)-1)*ny+j2(:);
        iH(index + (1:numel(e2))) = e1;
        jH(index + (1:numel(e2))) = e2;
        sH(index + (1:numel(e2))) = max(0,rmin-sqrt((i1-i2(:)).^2+(j1-j2(:)).^2));
        index = index + numel(e2);
    end
end
H = sparse(iH,jH,sH); Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,ny,nx);
xPhys = x;
iter = 0;
change = 1;
% micro
HoleRadius = 9;
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
%% make dir
str1 = strcat('Micro_Results_',num2str(EleNumPerRow),'x',num2str(EleNumPerCol),'_', num2str(vol_mi));
parent_dir_name =strcat('Experiment_results_', num2str(vol_ma),'\');
Micro_result_dir = strcat(parent_dir_name,str1);
mkdir(Micro_result_dir);
str2 = strcat('Macrostructure',num2str(nx),'x',num2str(ny),'_v',num2str(vol_ma));
Macro_results_dir = strcat(parent_dir_name,str2);
mkdir(Macro_results_dir);

%% START ITERATION
while change > 0.01
  iter = iter + 1;
  %% homogenization
  [dCH11,dCH12,dCH13,dCH22,dCH23,dCH33,Ke,A0,CH] = FEACalculation_CM(Phi,nelx,nely,EW,EH,EW,EH,E0,nu);
  %% FE-ANALYSIS
  sK = reshape(Ke(:)*max(xPhys(:),1e-3)'.^penal,64*nx*ny,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*Ke).*U(edofMat),2),ny,nx);
  c = sum(sum(x.^penal.*ce));
  dc = -penal*xPhys.^(penal-1).*ce;
  dv = ones(ny,nx);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  %% 
  [x,xPhys,change] = OC(x,dc,dv,nx,ny,vol_ma,ft);
  %%
  [~,Phi,Coef] = TopMMALSM(EleNumPerRow, EleNumPerCol,CentXYi,iter, DomainSize,vol_mi, passive,...
    nelx, nely, Coef, Phi, A, rmin, aspect,cIMQ,A0,reshape(Ke(:),8,8),U,nx,ny,dCH11,dCH12,dCH13,dCH22,dCH23,dCH33,EW,EH,CH,edofMat,xPhys,Micro_result_dir);

  
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',iter,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  figure(11);colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
end
%% OC update
function [x,xPhys,change] = OC(x,dc,dv,nx,ny,volfrac,ft)
l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(1e-3,max(x-move,min(1,min(x+move,x.*sqrt(max(1e-10,-dc./dv/lmid))))));
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
    end
    if sum(xPhys(:)) > volfrac*nx*ny, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

