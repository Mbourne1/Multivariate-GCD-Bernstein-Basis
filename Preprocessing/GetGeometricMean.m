function [GM_fx] = GetGeometricMean(fxy, n1_k1, n2_k2)
% Given the polynomial f(x,y) and g(x,y), get the geometric mean of their
% entries in the Sylvester matrix S_{k_{1},k_{2})(f,g).
%
% Inputs
%
% fxy_mtrx : (Matrix) Coefficients of the polynomial f(x,y)
%
% n1_k1 : (Int) Subscript of matrix T_{n1-k1,n2-k2}, the first partition of S_{k1,k2}
%
% n2_k2 : (Int) Subscript of matrix T_{n1-k1,n2-k2}, the first partition of S_{k1,k2}
%
% Outputs
%
% GM_fx : (Float) Geometric mean of entries in T_{n1-k1,n2-k2}(f)

% Build the matrix D^{-1}T_{n_{1}-k_{1},n_{2}-k_{2}}Q, the first partition
% of the Sylvester subresultant matrix
DT1Q1 = BuildSubresultant_Partition_Bivariate_2Polys(fxy, n1_k1, n2_k2);

% Get geometric mean of the nonzero entries
v1 = reshape(abs(DT1Q1(DT1Q1~=0)),1,[]);

GM_fx = geomean(v1);


end