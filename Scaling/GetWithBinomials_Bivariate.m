function fxy_bi = GetWithBinomials_Bivariate(fxy)
% Given the matrix of coefficients f(x,y), include the binomial
% coefficients in the matrix.
%
% Inputs
%
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% Outputs.
%
%
% fxy_bi : (Matrix) Coefficients of the polynomial f(x,y) including
%          binomial coefficients, so fxy_bi is the polynomial coefficients in the
%          modified Bernstein Basis.

% Get the degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get binomial coefficients.
bi_m1 = GetBinomials(m1);
bi_m2 = GetBinomials(m2);

% Get f(x,y) with binomials included (in modified Bernstein Basis).
fxy_bi = diag(bi_m1) * fxy * diag(bi_m2);

end
