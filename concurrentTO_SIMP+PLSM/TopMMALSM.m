function [fval, Phi0, coeff] = TopMMALSM(EleNumPerRow, EleNumPerCol,CentXYi,itte, DomainSize,volfrac,passive, ...
              nelx, nely, coeff, Phi0, A, rmin, aspect,cIMQ,A0,KE,UI,nx,ny,dCH11,dCH12,dCH13,dCH22,dCH23,dCH33,EW,EH,CH,edofMat,x,Micro_result_dir,nc)
% Start of MAIN program
% volfrac = vol;   % volume constraint
x = reshape (x,ny,nx);
tao = [1e-2];
tip =0;

format short g;
xval = coeff;    % Initialization of design variables:  (nely+1)*(nelx+1)-by-1-by-1

    % Sensitivity Analysis
    [df0dx,df0dx2, fval, dfdx, dfdx2, tfem, trbf] = IniParaMMA_CM(CentXYi, nelx, nely,...
               DomainSize, coeff, Phi0, A, cIMQ,volfrac, rmin, aspect,A0,KE,UI,nx,ny,dCH11,dCH12,dCH13,dCH22,dCH23,dCH33,EW,EH,edofMat,x);

% STEP 2 --  Finite Element Calculation and Shape Sensitivity Analysis
%
Tfem=[];Trbf=[];Tit=[];
%

% % Archive location of output files
warning off;
% str = strcat('Microstructure',num2str(EleNumPerRow),'x',num2str(EleNumPerCol));
% % str = strcat(str,'_v',num2str(volfrac));
% parent_dir_name =strcat('Experiment_results_',num2str(nx),'x',num2str(ny),'\');
% Micro_results_dir = strcat(parent_dir_name,str);
% mkdir(Micro_results_dir);
% str = strcat(Micro_result_dir,'\');
% str1 = strcat('Micro_cat_',num2str(nc));
% Micro_cat = strcat(str,str1);
% mkdir(Micro_cat);

        lower = 2*(min(coeff));
        upper = 2*(max(coeff));
        %
        [xmma]=OC(xval, A, lower, upper, volfrac, nelx, nely, df0dx, dfdx,itte,passive);
        
        % STEP 4 -- MMA optimization and Level set surface renewal
        %
        %[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = GCMMA(m,n,IterNum,xval,xmin,...
        %              xmax,xold1,xold2,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);        
        % xmma = xval - 1e-6 * df0dx;
        % Iteration results for next new iteration
        xval = xmma;       % (nely+1)*(nelx+1)-by-1
        coeff = xval;         % (nely+1)*(nelx+1)-by-1
        % "A" is the matrix comprosed by RBFs
        [R1,R2]=dwt(coeff,'db1');
        coeffx=[R1;R2];
        Dx=A*coeffx;
        f=idwt(Dx(1:length(R1)),Dx((length(R1)+1):end),'db1');
        f=f(1:(nelx+1)*(nely+1)); 
        Phi0 = f;               % (nely+1)*(nelx+1)-by-1             
        Phi0 = reshape(Phi0,nely+1,nelx+1);
        Phi0(find(passive)) = 1/3*max(max(Phi0));
        Phi0 = reshape(Phi0,(nely+1)*(nelx+1),1);

   % STEP 5 -- Results visualization to print optimal results.
    ft1 = 0;
    
    f0 = volfrac+fval;
    cMQ = cIMQ; 

    
    VisualResult(Phi0, f0,DomainSize,...
         itte, ft1, tfem, trbf,Tfem,Trbf, Tit,df0dx,xval,cMQ,...
         tao,EleNumPerCol,EleNumPerRow,Micro_result_dir,CH,itte);
      
    tip = tip+1;    
 
    % for restard calculation
    Resultsall(1, tip)= itte;
    
    optiresObj(1, tip)= itte;
    optiresVol(1, tip)= itte;
    optiresVol(2, tip)= fval(1);
    
   
    

% The end of the MAIN program  
    
% The end of the first-level sub-function "TopMMALSM(...)" 
