function [lambda,mu] = GetGeometricMean(fxy_mtrx,gxy_mtrx,k1,k2)
% Given teh polynomial f(x,y) and g(x,y), get the geometric mean of their
% entries in the Sylvester matrix S_{k_{1},k_{2})(f,g).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                     Inputs
%
%
% fxy_mtrx : Coefficients of the polynomial f(x,y)
%
% gxy_mtrx : Coefficients of the polynomial g(x,y)
%
% k1 : Degree of common divisor d(x,y)
%
% k2 : Degree of common divisor d(x,y)
%
%
%
%                     Outputs
%
% lambda : Geometric mean of entries in T(f)
%
% mu : Geometric mean of entries in T(g)


% Get degree of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_mtrx);

% Get number of columns in first partition
colsA = (n1-k1+1) * (n2-k2+1);

% Get number of columns in second partition
%colsB = (m1-k1+1) * (m2-k2+1);

Sk = BuildSubresultant(fxy_mtrx,gxy_mtrx,k1,k2,1,1,1);

C_f = Sk(:,1:colsA);
C_g = Sk(:,colsA+1:end);


lambda = geomean(reshape(abs(C_f(C_f~=0)),1,[]));
mu = geomean(reshape(abs(C_g(C_g~=0)),1,[]));

end