function VisualResult(Phi, f0, DomainSize,...
    IterNum, ft1, tfem, trbf,Tfem,Trbf, Tit,df0dx,xval,cMQ,...
    tao,EleNumPerCol,EleNumPerRow,parent_dir_name,CH,itte)

 % Output window configuration 
 nely = EleNumPerCol;
 nelx = EleNumPerRow;
 EleWidth = DomainSize(1) / EleNumPerRow;
 EleHight = DomainSize(2) / EleNumPerCol;
 
 bdwidth = 5;
 topbdwidth = 30;
 set(0,'Units','pixels');
 scnsize = get(0,'ScreenSize');

 FigNum = IterNum;
 Phi = reshape(Phi(:), (nely+1), (nelx+1)); 
 

 thickPara = 1;
 aspect=DomainSize(2)/nely;
 xAxis = [0:aspect/thickPara:DomainSize(1)];
 yAxis = [0:aspect/thickPara:DomainSize(2)];
 coeff_map = reshape(xval,nely*thickPara+1,nelx*thickPara+1);
 df0dx_map = reshape(-df0dx,nely*thickPara+1,nelx*thickPara+1);
 
 % Display the optimum results

 
 
 % First figure
 h1 = figure(1);
 set(h1,'visible','on');
 set( gcf, 'Color', [1, 1, 1] );
% if max(Phi(:))*min(Phi(:)) < -1e-19    
     contourf([0:nelx], [0:nely], Phi, [0 0]);  % draw the contour
     colormap([0 0 1]);
 %end
 axis equal; axis off; hold off;
%  FileName = [parent_dir_name,'/contour_', int2str(FigNum),'.fig'];
%  saveas(h1,FileName);
FileName = [parent_dir_name,'/contour_', int2str(FigNum),'.tif'];
 saveas(h1,FileName);
 
%  h2 = figure(2);
%  set(h2,'visible','on');
%  set( gcf, 'Color', [1, 1, 1] );
% % if max(Phi(:))*min(Phi(:)) < -1e-19    
%      contourf([0:nelx], [0:nely], Phi, [0 0]);  % draw the contour
%      colormap([0,0,0]);
%  %end
%  axis equal; axis off; hold off;
%  FileName = [parent_dir_name1,'/contour_', int2str(FigNum),'.fig'];
%  saveas(h2,FileName);
% FileName = [parent_dir_name1,'/contour_', int2str(FigNum),'.tif'];
%  saveas(h2,FileName);
 
  % Third figure
 % h3 = figure(3);
 % set(h3,'visible','on');
 % set( gcf, 'Color', [1, 1, 1] );
 % % if max(Phi(:))*min(Phi(:)) < -1e-19    
 %     contourf([0:nelx], [0:nely], Phi, [0 0]);  % draw the contour
 % %end
 % alpha(0.2);
 % hold on;
 % h3_surf = surf([0:nelx], [0:nely], Phi);
 % view([165,120]);
 % set(h3_surf,'FaceLighting','phong','FaceColor','interp',...
 %     'AmbientStrength',1.0,'DiffuseStrength',1);
 % light('Position', [1 1 1], 'Style','infinite');
 % axis off;
 % hold off; %axis equal; 
%  FileName = [parent_dir_name,'/surf_',int2str(FigNum),'.fig'];
% %  saveas(h3, FileName);
%  FileName = [parent_dir_name,'/surf_',int2str(FigNum),'.tif'];
%  saveas(h3, FileName);
%  % 
 ft2 = cputime - ft1;
 Tfem = [Tfem; tfem];
 Trbf = [Trbf; trbf];
 Tit = [Tit; ft2];
  
 DH=reshape(CH,3,3);
 FileName=[parent_dir_name,'\DH.txt'];
 save(FileName, 'DH','-ascii');
 
  FileName=[parent_dir_name,'\Micro_Phi', '_',int2str(itte)];
 save(FileName, 'Phi');
 
end
