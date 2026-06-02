function [deltaPhi] = deltaWholePhi(Phi, aspect)
[x,y] = size(Phi);
epsilon = 2*aspect;
for i = 1:x;
    for j = 1:y;
       deltaPhi(i,j) = Phi(i,j)^2/(epsilon^2+Phi(i,j)^2)/pi;
       %deltaPhi(i,j) = epsilon/(epsilon^2+Phi(i,j)^2)/pi;
    end
end