% Start the "OC" subfunction
function [coeffnew]=OC(coeff, A, lower, upper, volfrac, nelx, nely, df0dx, dfdx,itte,passive) 

% nely = 50;nelx=100;
% coeff = 0.6*ones((nely+1)*(nelx+1),1);
% lower = 0.1;
% upper = 1;
% df0dx = 0.025*ones((nely+1)*(nelx+1),1);
% dfdx = 0.275*ones((nely+1)*(nelx+1),1);

xold = coeff;
move = 0.005; damp = 0.5;
% if itte < 12
%     move = 0.01;damp = 0.5;
% else 
%     move = 0.06;damp = 0.40;
% end
% Start -- Bi-section method to update the Lagrange Multiplier
l1 = 0; 
l2 = 100000;
xold = reshape(xold, nely+1, nelx+1);
coeff = reshape(coeff, nely+1, nelx+1);
df0dx = reshape(df0dx, nely+1, nelx+1);
dfdx = reshape(dfdx, nely+1, nelx+1);
lmid = 0.5*(l2+l1);
while (l2-l1)/(l2+l1) > 1e-5 && l2 > 1e-40    
    % Converting the design variables belonging to [0 1]
    x = (xold - lower)./(upper - lower); 
    xnew = max(0.001, max(x - move, min(1., min(x + move,...
                x.*(max(1e-12, -df0dx./(max(1e-12,dfdx.*lmid))).^damp)))));  
    % Converting the design variables belonging to [lower  upper]
    x1 = (xnew.*(upper - lower) + lower);    
    coeffnew1 = reshape(x1, (nely+1)*(nelx+1), 1); 
    [R1,R2]=dwt(coeffnew1,'db1');
    coeffnew=[R1;R2];
    tempphi = A*coeffnew;
    tempphi0=idwt(tempphi(1:length(R1)),tempphi((length(R1)+1):end),'db1');
    tempphi=tempphi0(1:(nelx+1)*(nely+1));coeffnew=coeffnew1;
        TT = reshape(tempphi,nely+1,nelx+1);
        TT(find(passive)) = 1/3*max(max(TT));
        tempphi = reshape(TT,(nely+1)*(nelx+1),1);
    [A0, Arc, den, xy0] = BoundEle(reshape(tempphi, nely+1, nelx+1)); 
    if A0 - volfrac*(nely*nelx) > 0 %|| A0 - 0.2*(nely*nelx) < 0
        l1 = lmid;
    else
        l2 = lmid;
    end
    lmid = 0.5*(l2+l1);
end
% End the "OC" subfunction
