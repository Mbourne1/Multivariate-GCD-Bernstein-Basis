function [lambda] = GetGeometricMean(fxy_matrix,n1_k1,n2_k2)
% Given teh polynomial f(x,y) and g(x,y), get the geometric mean of their
% entries in the Sylvester matrix S_{k_{1},k_{2})(f,g).
%
% Inputs
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
%                     Outputs
%
%
% lambda : Geometric mean of entries in T(f)



% Get number of columns in second partition
%colsB = (m1-k1+1) * (m2-k2+1);
DT1Q1 = BuildDT1Q1(fxy_matrix,n1_k1,n2_k2);

lambda = geomean(reshape(abs(DT1Q1(DT1Q1~=0)),1,[]));


end