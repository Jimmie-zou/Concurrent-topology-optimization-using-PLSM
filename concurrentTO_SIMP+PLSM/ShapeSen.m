function [dc] = ShapeSen(dCH11,dCH12,dCH13,dCH22,dCH23,dCH33,KE,UI,nx,ny,NN,EW,EH,edofMat,x)

c = 0;
dc = zeros(1,NN);
for i = 1:NN
    dD = [dCH11(i) dCH12(i) dCH13(i); dCH12(i) dCH22(i) dCH23(i); dCH13(i) dCH23(i) dCH33(i)];
    [dKE] = StiffnessMatrix(dD,EW,EH);
       
    vaa = reshape(sum((UI(edofMat)*dKE).*UI(edofMat),2),ny,nx);
    va=sum(sum(x.*vaa));
    dc(i) = va;
end
% ce = reshape(sum((UI(edofMat) * KE) .* UI(edofMat), 2), ny, nx);
% c = sum(ce, 'all');
