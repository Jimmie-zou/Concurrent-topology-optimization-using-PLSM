function [df0dx,df0dx2, fval, dfdx, dfdx2, tfem, trbf] = IniParaMMA_CM(CentXYi, nelx, nely,...
    DomainSize, coeff, Phi, A, cIMQ,volfrac, rmin, aspect,A0,KE,UI,nx,ny,dCH11,dCH12,dCH13,dCH22,dCH23,dCH33,EW,EH,edofMat,x_core)


f0val   = 0;  
df0dx   =coeff*0;
dfdx    =coeff*0;
df0dx2  =coeff*0;
dfdx2   =coeff*0;

% Call "deltaWholePhi"
deltaPhi = deltaWholePhi(Phi,aspect);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Start the sub-function "IniParaMMA_CM(...)"
tao = 1e-2;
c = cIMQ;
t1 = cputime; 
gridcoor = CentXYi;
[nKnots  dim] = size(Phi);
centnum = nKnots;
EleNumPerRow = nelx;
EleNumPerCol = nely;
EleWidth = DomainSize(1) / EleNumPerRow;
EleHight = DomainSize(2) / EleNumPerCol;

tfem = cputime-t1;

% CALL sub-function to perform FEA analysis
%[CH,dc,c, xy0, Arc, A0, den] = FEACalculation_CM1(x0, EleWidth, EleHight, nely, nelx);
NN = (nely+1)*(nelx+1);
[dc] = ShapeSen(dCH11,dCH12,dCH13,dCH22,dCH23,dCH33,KE,UI,nx,ny,NN,EW,EH,edofMat,x_core);

dSE_dAlpha = zeros((nelx+1)*(nely+1), 1);
dVolume_dAlpha = zeros((nelx+1)*(nely+1), 1);

ell = 0;
dc = reshape(dc,nely+1,nelx+1);
dc=(rot90(dc,2)+dc)/2;

dc(:,1:(nelx+1+1)/2-1) = (dc(:,end:-1:(nelx+1+1)/2+1) + dc(:,1:(nelx+1+1)/2-1))/2;
dc(:,end:-1:(nelx+1+1)/2+1) = dc(:,1:(nelx+1+1)/2-1);
dc(1:(nely+1+1)/2-1,:) = (dc(end:-1:(nely+1+1)/2+1,:) + dc(1:(nely+1+1)/2-1,:))/2;
dc(end:-1:(nely+1+1)/2+1,:) = dc(1:(nely+1+1)/2-1,:);

dc = reshape(dc,1,(nelx+1)*(nely+1));
eSE = dc;
[eSE] = Check(nelx+1,nely+1,rmin,eSE);     
% Bi-section to calculate the Lagrange multiplier if necessary

% End the first step -- formulating the objective function
%
trbf=cputime-t1;
del = [deltaPhi(:);0];
[R1,R2] = dwt(del,'db1');
deltaPhix = [R1;R2];
DeltaPhix = A*deltaPhix;
DeltaPhi = idwt(DeltaPhix(1:length(R1)),DeltaPhix((length(R1)+1):end),'db1');
DeltaPhi = DeltaPhi(1:(nelx+1)*(nely+1));
%
% Objective function and design sensitivity.
%
dObj_dAlpha = -DeltaPhi(:).*eSE(:).*aspect*aspect;

df0dx = dObj_dAlpha;
df0dx2 = 0*df0dx;
%
% Constraint and design sensitivity.
%
vol = A0 - volfrac; 
dVol_dAlpha = DeltaPhi(:).*aspect*aspect;

fval = [vol]';
dfdx = [dVol_dAlpha]';
dfdx2 = 0*dfdx;
% End -- calculating the design functions and their sensitivities
