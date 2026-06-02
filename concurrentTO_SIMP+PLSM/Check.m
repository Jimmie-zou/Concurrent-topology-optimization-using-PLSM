function [dcn] = Check(nelx,nely,rmin,dc)
% Filter for Sensitivity
% dcn=zeros(nely+1,nelx+1);
% den=zeros((nely+1)*(nelx+1),2);
% for ely = 1:nely
%     for elx = 1:nelx
%         n1 = (nely+1)*(elx-1)+ely;
%         n2 = (nely+1)* elx   +ely;
%         den(n1,1)=x(ely,elx)+den(n1,1);
%         den(n1,2)=1+den(n1,2);
%         den(n1+1,1)=x(ely,elx)+den(n1+1,1);
%         den(n1+1,2)=1+den(n1+1,2);
%         den(n2,1)=x(ely,elx)+den(n2,1);
%         den(n2,2)=1+den(n2,2);
%         den(n2+1,1)=x(ely,elx)+den(n2+1,1);
%         den(n2+1,2)=1+den(n2+1,2);
%     end
% end
% den(:,1)=den(:,1)./den(:,2);
% den=reshape(den(:,1),nely+1,nelx+1);

% for i = 1:nelx
%     for j = 1:nely
%         sum=0.0;
%         for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
%             for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
%                 fac = rmin-sqrt((i-k)^2+(j-l)^2);
%                 sum = sum+max(0,fac);
%                 dcn(j,i) = dcn(j,i) + max(0,fac)*den(l,k)*dc(l,k);
%             end
%         end
%         dcn(j,i) = dcn(j,i)/(den(j,i)*sum);
%     end
% end

% Filter 
xcon = reshape(dc,nely,nelx);
xconm=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sumcf=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
          rr = rmin - sqrt((i-k)^2+(j-l)^2);
          rrr = 1 - sqrt((i-k)^2+(j-l)^2)/rmin;
          convf  = (3/(pi*rmin.^2))*max(0,rrr)*1.0;          
          sumcf = sumcf + max(0,convf);          
          xconm(j,i) = xconm(j,i) + max(0,convf)*xcon(l,k); 
      end
    end 
    xconm(j,i) = xconm(j,i)/(sumcf);
  end
end
dcn = reshape(xconm,1,nely*nelx);
% The end of this function