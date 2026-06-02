function [Ke] = StiffnessMatrix(D, EW,EH)
% D - Elastic material matrix ;
Ke = zeros(8, 8); 
% Loop 4 Gauss integration points
s = [ -0.5773;  0.5773]/EW; 
t = [ -0.5773;  0.5773]/EH;  
h = 1;
W = [1 1; 1 1]; 
for i = 1:2 
    for j = 1:2
        % Calculate strain matrix for one Guass point
        [Bu, J]=StrainMarix(h,s(i),t(j),EW,EH);   
        Ke = Ke + (Bu'*D*Bu)*det(J)*h*W(i,j)*W(i,j); % 8*8
%         Me = Me + mass_density*(N'*N)*det(J)*h*W(i,j)*W(i,j);
    end
end


function [Bu, J] = StrainMarix(h, s, t,EW,EH)
R = 1/4.*[ (-1+t) (1-t) (1+t) (-1-t) ; (-1+s) (-1-s) (1+s) (1-s) ];
% Jacobian of 2D transform connecting (x,y) with (xi,niu) 
Coords = [ 0   0; EW  0 ;EW  EH; 0  EH];
J = R*Coords;
dN =inv(J)*R;
% Small linear displacement strain-displace transfer matrix
B1 = [dN(1,1)   0         dN(1,2)   0         dN(1,3)    0              dN(1,4)    0           ];
B2 = [0         dN(2,1)   0         dN(2,2)   0               dN(2,3)   0               dN(2,4)];
B3 = [dN(2,1)   dN(1,1)   dN(2,2)   dN(1,2)   dN(2,3)     dN(1,3)   dN(2,4)    dN(1,4)];
Bu = [B1;   B2;   B3]; 
%--------------------------------------------------------------------------
