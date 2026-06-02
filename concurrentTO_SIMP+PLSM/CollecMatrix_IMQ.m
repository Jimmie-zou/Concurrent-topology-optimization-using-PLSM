function [A] = CollecMatrix_IMQ(NodeXYi, CentXYi, cIMQ, DomainSize, nelx, nely)
% to construct  A -- (n+dim+1)*(n+dim+1) matrix
c = cIMQ;
[n dim] = size(CentXYi);   % n-by-2; n = (nely+1)(nelx+1), dim = 2.
[m dim] = size(NodeXYi); % n-by-2; n = (nely+1)(nelx+1), dim = 2.
% Method-1. IMQRBF
A = zeros(n, n);
for i = 1:n
    for j = 1:i
        r = norm(CentXYi(i,:) - CentXYi(j,:));
        temp = max(exp(-r*r*c),1e-320);
        %temp = sqrt(r*r - c*c);% Inverse multiquadric RBF (IMQ)
        A(i,j) = temp;
        A(j,i) = temp; % A(j,i) means the value of the basis function of knot i at point x_j -- the jth grid node. 
    end
end

Cx = A;
B=[Cx,zeros((nelx+1)*(nely+1),1)];
B=[B;zeros(1,(nelx+1)*(nely+1)+1)];
[cc,s]=wavedec2(B,1,'db1');
thr=sum(abs(cc))/(s(3,2)*s(3,1));
cc(abs(cc)<(1)*thr)=0;cc=sparse(cc);
[CA]=reshape(cc(1:s(1,1)*s(1,1)),s(1,1),s(1,1));
[CH]=reshape(cc(s(1,1)*s(1,1)+1:2*s(2,1)*s(2,1)),s(2,1),s(2,1));
[CV]=reshape(cc(2*s(2,1)*s(2,1)+1:3*s(1,2)*s(1,2)),s(1,2),s(1,2));
[CD]=reshape(cc(3*s(1,2)*s(1,2)+1:4*s(2,2)*s(2,2)),s(2,2),s(2,2));
A=[CA,CV;CH,CD];A(end)=0;
