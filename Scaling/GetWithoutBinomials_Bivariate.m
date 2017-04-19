function [fxy] = GetWithoutBinomials_Bivariate(fxy_bi)
% Given coefficients of the polynomial f(x,y) in scaled Bernstein basis,
% return the coefficients of the polynomial in the Bernstein basis.
%
% % Inputs
%
% fxy_bi : (Matrix) Coefficients of f(x,y) in scaled Bernstein Basis.
%
% % Outputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) in Bernstein basis


% Get degree structure of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy_bi);

% Get the binomials of basis elements B_{i}^{m_{1}}(x) ie
% \nchoosek(m_{1},i)
bi_m1 = GetBinomials(m1);
mat1 = diag(1./bi_m1);

% Get the binomials of basis elements B_{i}^{m_{2}}(y) ie
% \nchoosek(m_{2},i)
bi_m2 = GetBinomials(m2);
mat2 = diag(1./bi_m2);

% Multiply coefficients of f(x,y) by binomials
fxy = mat1 * fxy_bi * mat2;

end