function [lambda] = GetGeometricMean(fxy_matrix,n1_k1,n2_k2)
% Given the polynomial f(x,y) and g(x,y), get the geometric mean of their
% entries in the Sylvester matrix S_{k_{1},k_{2})(f,g).
%
% Inputs
%
% fxy_mtrx : Coefficients of the polynomial f(x,y)
%
% n1_k1 : Degree of common divisor v(x,y)
%
% n2_k2 : Degree of common divisor v(x,y)
%
% Outputs
%
% lambda : Geometric mean of entries in T(f)

% Build the matrix D^{-1}T_{n_{1}-k_{1},n_{2}-k_{2}}Q, the first partition
% of the Sylvester subresultant matrix
DT1Q1 = BuildDT1Q1(fxy_matrix,n1_k1,n2_k2);

% Get geometric mean of the nonzero entries
lambda = geomean(reshape(abs(DT1Q1(DT1Q1~=0)),1,[]));


end