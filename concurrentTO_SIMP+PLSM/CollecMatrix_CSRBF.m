function [A] = CollecMatrix_CSRBF( CentXYi, DomainSize, nelx, nely)

dmax = 3.5; %% CSBRF插值的半紧支域径，即单元i周边距离3.5以内的单元才有效 %%%
numgrids = (nelx+1)*(nely+1);
xspace = DomainSize(1)/nelx;
yspace = DomainSize(2)/nely;
%% method 1
dm(1:numgrids,1) = dmax*xspace*ones(numgrids,1);
dm(1:numgrids,2) = dmax*yspace*ones(numgrids,1);
% A = sparse(numgrids, numgrids);
% gridcoor = CentXYi;
% for i = 1:numgrids
%     dmI = sqrt(dm(i, 1).^2 + dm(i, 2).^2);
%     for j = 1:numgrids
%         dI = norm(gridcoor(i, :) - gridcoor(j, :)); %% 节点i到设计域内所有节点的距离 %%%
%         r  = dI./dmI; %% 在紧支域半径内的，距离r小于1，通过紧支径向基函数插值可得到有效的值，%%%
%         %%%% 在紧支域半径外的，距离r大于1，通过紧支径向基函数插值后为0 %%%
%         csrbfrC2 = (max(0, (1-r)).^4).*(4*r + 1);
% %         %csrbfrC4 = (max(0, (1-r)).^6).*(35*r.^2 + 18*r + 3);
% %         %csrbfrC6 = (max(0, (1-r)).^8).*(32*r.^3 + 25*r.^2 + 8*r + 1);
%         A(j, i) = csrbfrC2;  %%% 紧支径向基函数值，A维度为 (nelx+1)*(nely+1) by (nelx+1)*(nely+1)
%         %%%% 每一列即为一个节点到设计域内所有节点的距离，经CSRBF插值后的值，总共(nelx+1)*(nely+1)列 %%%
%     end
% end

%% method 2
dmI=sqrt(dm(:,1).^2+dm(:,2).^2);
grid_x_0 = CentXYi(:,1);
grid_y_0 = CentXYi(:,2);
grid_x_1 = repmat(grid_x_0,1,numgrids);
grid_x_2 = repmat(grid_x_0',numgrids,1);
grid_y_1 = repmat(grid_y_0,1,numgrids);
grid_y_2 = repmat(grid_y_0',numgrids,1);
grid_xy = sqrt((grid_x_1-grid_x_2).^2 + (grid_y_1-grid_y_2).^2);
r = grid_xy./dmI;
A0 = (max(0, (1-r)).^4).*(4*r + 1);
A = sparse(A0);
