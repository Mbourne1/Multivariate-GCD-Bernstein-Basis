function [GM_f,GM_g] = getGeometricMean(fxy_mtrx,gxy_mtrx,k1,k2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                     Inputs


% fxy_mtrx_bi

% gxy_mtrx_bi

% k1

% k2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                     Outputs




[rows,cols] = size(fxy_mtrx);
m1 = rows -1;
m2 = cols -1;

[rows,cols] = size(gxy_mtrx);
n1 = rows -1;
n2 = cols -1;

colsA = (n1-k1+1) * (n2-k2+1);
colsB = (m1-k1+1) * (m2-k2+1);

Sk = BuildSubresultant(fxy_mtrx,gxy_mtrx,k1,k2,1,1);

C_f = Sk(:,1:colsA);
C_g = Sk(:,colsA+1:end);



GM_f = geomean(reshape(abs(C_f(C_f~=0)),1,[]));
GM_g = geomean(reshape(abs(C_g(C_g~=0)),1,[]));
end