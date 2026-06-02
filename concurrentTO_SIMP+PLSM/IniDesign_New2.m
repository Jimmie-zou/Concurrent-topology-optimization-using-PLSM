function [Phi0, RBFGridXYi, FEANodeXYi] = IniDesign_New2(nelx, nely, DomainSize, holeRadius);
% % 
% close all;
% clear all;
% nelx=60;nely=30;
% DomainSize=[60  30];
% %
EleNumPerRow=nelx;
EleNumPerCol=nely;
EleWidth = DomainSize(1) / EleNumPerRow;
EleHight = DomainSize(2) / EleNumPerCol;

HoleRadius = holeRadius*(DomainSize(1) / EleNumPerRow);
% With element density increasing to be slightly larger


% RBF grids
Gridti1 = 0 : EleWidth : DomainSize(1);
Gridti2 = 0 : EleHight : DomainSize(2);
[XI,YI] = meshgrid(Gridti1, Gridti2);
x = XI(:);
y = YI(:);
RBFGridXYi = [x   y];

% FEA nodes
Node_Xi = 0 : EleWidth : DomainSize(1);
Node_Yi = 0 : EleHight : DomainSize(2);
[NodeXi, NodeYi] = meshgrid(Node_Xi , Node_Yi);
FEANodeXi = NodeXi(:); FEANodeYi = NodeYi(:); 
FEANodeXYi = [FEANodeXi   FEANodeYi];



% cx = DomainSize(1)*[1/4  3/4  1/4  3/4   0/4  2/4  4/4 ];
% cy = DomainSize(2)*[0/2  0/2   1     1   1/2   1/2    1/2];
% cx = [cx(:)];
% cy = [cy(:)];
% cx = DomainSize(1)*[0 0 1 1 1/2];
% cy = DomainSize(2)*[0 1 0 1 1/2];
% cx = [cx(:)];
% cy = [cy(:)];
cx = DomainSize(1)*[0   1   1/2 1/2];
cy = DomainSize(2)*[1/2 1/2 0   1];
cx = [cx(:)];
cy = [cy(:)];



% cx = DomainSize(1)*[1/4; 1/4; 3/4];
% cy = DomainSize(2)*[1/4; 3/4; 2/4];
% cx = cx(:);
% cy = cy(:);

% cx = DomainSize(1)*[1/4  3/4];
% cy = DomainSize(2)*[1/2  1/2];
% cx = [cx(:)];
% cy = [cy(:)];

% xeven = DomainSize(1)*[1/4 : 2/4: 1];
% yeven = DomainSize(2)*[1/4 : 2/4: 1];
% xodd = DomainSize(1)*[0/2 : 2/4: 1];
% yodd = DomainSize(2)*[0/2 : 2/4: 1];

% [cx1, cy1] = meshgrid(xeven, yodd);
% [cx2, cy2] = meshgrid(xodd, yeven);
% cx = [cx1(:);cx2(:)];
% cy = [cy1(:);cy2(:)];

% Considering boundary 
for m = 1 : length( x )
    for k = 1 : length( cx )
        tmpz( k ) = -sqrt ( ( x(m) - cx ( k ) ) ^2 + ( y ( m ) - cy ( k ) ) ^2 ) + HoleRadius;
        tmpz( k+1) = -x(m);
        tmpz( k+2) = -y(m);
        tmpz( k+3) = -(DomainSize(1)-x(m));
        tmpz( k+4) = -(DomainSize(2)-y(m));
    end; 
    Phi(m)=max(tmpz);
end;
Phi = - Phi.';
x0 = reshape(Phi,nely+1,nelx+1);
for i=1:(nelx+1)
    for j=1:(nely+1)
        xx(nely+1-j+1,i)=x0(j,i);
    end
end
Phi0=reshape(xx,(nelx+1)*(nely+1),1);
% % 
% % % % 
% Phitest=reshape(Phi0,(nely+1),(nelx+1));
% last = reshape(Phi, nely+1, nelx+1);
% 
% h1 = figure(1);
% set( gcf, 'Color', [1, 1, 1] );
% contourf([0:nelx], [0:nely], last, [0 0]); 
% colormap([0.1 0.1 0.5]); 
% axis off;
% axis equal;
% hold on;
% 
% h2 = figure(2);
% set( gcf, 'Color', [1, 1, 1] );
% h_surf=surf([0:nelx], [0:nely], last);
% view([165,120]);
% set(h_surf,'FaceLighting','phong','FaceColor','interp',...
%      'AmbientStrength',1.0,'DiffuseStrength',1);
% light('Position', [1 1 1], 'Style','infinite');
% axis off;
% hold off;
% % The end of this function.