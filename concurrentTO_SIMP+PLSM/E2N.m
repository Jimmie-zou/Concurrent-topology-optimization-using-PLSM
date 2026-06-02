function [N] = E2N(nelx, nely, E)
% To circulate the "eSE" to all knots
% nely = 2; nelx = 4; 
% E = 2.5*ones(nely*nelx,1);
E = reshape(E, nely, nelx);
N = zeros((nely+1)*(nelx+1),1);
n = zeros((nely+1)*(nelx+1),1);
for elx = 1:nelx
    for ely = 1:nely
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx   +ely;
        n0=[n1;n2;n2+1;n1+1];
        for i=1:4,
            N(n0(i),1) = N(n0(i),1) + E(ely,elx);
            n(n0(i),1)=n(n0(i),1)+1;
        end
    end
end
N = N./n;
N = reshape(N, nely+1, nelx+1);