function [A0, Arc, den, xy] = BoundEle(xx)
%                   /\  y
%                   |
%                   -------     
%                   |     |
%                   |     |
%                   ------------>  x
% this is the local coordinates in one element 
%  "xx" refers to matrix-> "Phi0(nely+1,nelx+1)"
% Solid Area A0

A0 = 0;
Arc = 0;
[ny, nx] = size(xx);
nelx = nx - 1;
nely = ny - 1;
% local coordinates of the boundary line
xy = zeros(nelx*nely, 4);
% element area (density)
den = zeros(nely, nelx);
% boundary nodes
for i=1:nely,
    for j=1:nelx,
        e1 = nely*(j -1)+i;%element number;
        % Small Triangle Case 1
        if xx(i,j)>0 & xx(i+1,j)<0 & xx(i,j+1)<0 & xx(i+1,j+1)<0
            a = -xx(i,j)/(xx(i,j+1)-xx(i,j));
            b = -xx(i,j)/(xx(i+1,j)-xx(i,j));
            Arc = Arc + sqrt(a^2+b^2);
            A0 = A0 + a*b/2.0;
            den(i, j) = a*b/2.0;
            xy(e1,1) = a; xy(e1,2) = 1;
            xy(e1,4) = 1-b;

        % Small Triangle Case 2
        elseif xx(i+1,j)>0 & xx(i,j)<0 & xx(i,j+1)<0 & xx(i+1,j+1)<0
            a=-xx(i+1,j)/(xx(i+1,j+1)-xx(i+1,j));
            b=-xx(i+1,j)/(xx(i,j)-xx(i+1,j));
            Arc=Arc+sqrt(a^2+b^2);
            A0=A0+a*b/2.0;
            den(i,j)=a*b/2.0;
            xy(e1,1)=a;
            xy(e1,4)=b;

        %Small Triangle Case 3
        elseif xx(i+1,j+1)>0 & xx(i+1,j)<0 & xx(i,j+1)<0 & xx(i,j)<0
            a=-xx(i+1,j+1)/(xx(i+1,j)-xx(i+1,j+1));
            b=-xx(i+1,j+1)/(xx(i,j+1)-xx(i+1,j+1));
            Arc=Arc+sqrt(a^2+b^2);
            A0=A0+a*b/2.0;
            den(i,j)=a*b/2.0;
            xy(e1,1)=1-a;
            xy(e1,4)=b;xy(e1,3)=1;

        %Small Triangle Case 4
        elseif xx(i,j+1)>0 & xx(i+1,j)<0 & xx(i,j)<0 & xx(i+1,j+1)<0
            a=-xx(i,j+1)/(xx(i,j)-xx(i,j+1));
            b=-xx(i,j+1)/(xx(i+1,j+1)-xx(i,j+1));
            Arc=Arc+sqrt(a^2+b^2);
            A0=A0+a*b/2.0;
            den(i,j)=a*b/2.0;
            xy(e1,1)=1-a;xy(e1,2)=1;
            xy(e1,4)=1-b;xy(e1,3)=1;

        %Triangle Case 1
        elseif xx(i,j)<0 & xx(i+1,j)>0 & xx(i,j+1)>0 & xx(i+1,j+1)>0
            a=-xx(i,j)/(xx(i,j+1)-xx(i,j));
            b=-xx(i,j)/(xx(i+1,j)-xx(i,j));
            Arc=Arc+sqrt(a^2+b^2);
            A0=A0+1-a*b/2.0;
            den(i,j)=1-a*b/2.0;
            xy(e1,1)=a;xy(e1,2)=1;
            xy(e1,4)=1-b;

        %Triangle Case 2
        elseif xx(i+1,j)<0 & xx(i,j)>0 & xx(i,j+1)>0 & xx(i+1,j+1)>0
            a=-xx(i+1,j)/(xx(i+1,j+1)-xx(i+1,j));
            b=-xx(i+1,j)/(xx(i,j)-xx(i+1,j));
            Arc=Arc+sqrt(a^2+b^2);
            A0=A0+1-a*b/2.0;
            den(i,j)=1-a*b/2.0;
            xy(e1,1)=a;
            xy(e1,4)=b;

        %Triangle Case 3
        elseif xx(i+1,j+1)<0 & xx(i+1,j)>0 & xx(i,j+1)>0 & xx(i,j)>0
            a=-xx(i+1,j+1)/(xx(i+1,j)-xx(i+1,j+1));
            b=-xx(i+1,j+1)/(xx(i,j+1)-xx(i+1,j+1));
            Arc=Arc+sqrt(a^2+b^2);
            A0=A0+1-a*b/2.0;
            den(i,j)=1-a*b/2.0;
            xy(e1,1)=1-a;
            xy(e1,4)=b;xy(e1,3)=1;

        %Triangle Case 4
        elseif xx(i,j+1)<0 & xx(i+1,j)>0 & xx(i,j)>0 & xx(i+1,j+1)>0
            a=-xx(i,j+1)/(xx(i,j)-xx(i,j+1));
            b=-xx(i,j+1)/(xx(i+1,j+1)-xx(i,j+1));
            Arc=Arc+sqrt(a^2+b^2);
            A0=A0+1-a*b/2.0;
            den(i,j)=1-a*b/2.0;
            xy(e1,1)=1-a;xy(e1,2)=1;
            xy(e1,4)=1-b;xy(e1,3)=1;

        %Quadrilateral Case 1
        elseif xx(i,j)<=0 & xx(i+1,j)<=0 & xx(i,j+1)>0 & xx(i+1,j+1)>0
            a=-xx(i,j)/(xx(i,j+1)-xx(i,j));
            b=-xx(i+1,j)/(xx(i+1,j+1)-xx(i+1,j));
            Arc=Arc+sqrt(a^2+b^2);
            A0=A0+1-(a+b)/2.0;
            den(i,j)=1-(a+b)/2.0;
            xy(e1,1)=a;xy(e1,2)=1;
            xy(e1,3)=b;

        %Quadrilateral Case 2
        elseif xx(i,j)>0 & xx(i+1,j)>0 & xx(i,j+1)<=0 & xx(i+1,j+1)<=0
            a=-xx(i,j+1)/(xx(i,j)-xx(i,j+1));
            b=-xx(i+1,j+1)/(xx(i+1,j)-xx(i+1,j+1));
            Arc=Arc+sqrt(a^2+b^2);
            A0=A0+1-(a+b)/2.0;
            den(i,j)=1-(a+b)/2.0;
            xy(e1,1)=1-a;xy(e1,2)=1;
            xy(e1,3)=1-b;

        %Quadrilateral Case 3
        elseif xx(i,j)<=0 & xx(i,j+1)<=0 & xx(i+1,j)>0 & xx(i+1,j+1)>0
            a=-xx(i,j)/(xx(i+1,j)-xx(i,j));
            b=-xx(i,j+1)/(xx(i+1,j+1)-xx(i,j+1));
            Arc=Arc+sqrt(a^2+b^2);
            A0=A0+1-(a+b)/2.0;
            den(i,j)=1-(a+b)/2.0;
            xy(e1,2)=1-a;
            xy(e1,4)=1-b;xy(e1,3)=1;

        %Quadrilateral Case 4
        elseif xx(i,j)>0 & xx(i,j+1)>0 & xx(i+1,j)<=0 & xx(i+1,j+1)<=0
            a=-xx(i+1,j)/(xx(i,j)-xx(i+1,j));
            b=-xx(i+1,j+1)/(xx(i,j+1)-xx(i+1,j+1));
            Arc=Arc+sqrt(a^2+b^2);
            A0=A0+1-(a+b)/2.0;
            den(i,j)=1-(a+b)/2.0;
            xy(e1,2)=a;
            xy(e1,4)=b;xy(e1,3)=1;

        %Diagonal Case 1 (Case 7)
        elseif xx(i,j)>0 & xx(i+1,j)==0 & xx(i,j+1)==0 & xx(i+1,j+1)<0
            Arc=Arc+sqrt(2);
            A0=A0+1/2.0;
            den(i,j)=1/2.0;
            xy(e1,3)=1;
            xy(e1,4)=1;

        %Diagonal Case 1 (Case -7)
        elseif xx(i,j)<0 & xx(i+1,j)==0 & xx(i,j+1)==0 & xx(i+1,j+1)>0
            Arc=Arc+sqrt(2);
            A0=A0+1/2.0;
            den(i,j)=1/2.0;
            xy(e1,3)=1;
            xy(e1,4)=1;

        %Diagonal Case 1 (Case 8)
        elseif xx(i,j)==0 & xx(i+1,j)>0 & xx(i,j+1)<0 & xx(i+1,j+1)==0
            Arc=Arc+sqrt(2);
            A0=A0+1/2.0;
            den(i,j)=1/2.0;
            xy(e1,2)=1;
            xy(e1,3)=1;

        %Diagonal Case 1 (Case -8)
        elseif xx(i,j)==0 & xx(i+1,j)<0 & xx(i,j+1)>0 & xx(i+1,j+1)==0
            Arc=Arc+sqrt(2);
            A0=A0+1/2.0;
            den(i,j)=1/2.0;
            xy(e1,2)=1;
            xy(e1,3)=1;
            
        % boundary line
        elseif i==1 | i==nely | j==1 | j==nelx,
            if xx(i,j)>=0 & xx(i+1,j)>=0 & xx(i,j+1)>=0 & xx(i+1,j+1)>=0
                Arc=Arc+1;
                A0=A0+1;
                den(i,j)=1;
            end

        %The whole elment solid
        elseif xx(i,j)>=0 & xx(i+1,j)>=0 & xx(i,j+1)>=0 & xx(i+1,j+1)>=0
            A0=A0+1;
            den(i,j)=1;
        end

    end
end


% end sub-function "(...)=gridBound(...)"