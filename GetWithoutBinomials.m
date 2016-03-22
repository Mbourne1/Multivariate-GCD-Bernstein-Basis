function [fxy] = GetWithoutBinomials(fxy_bi)

% Get degree structure of f(x,y)
[m1,m2] = GetDegree(fxy_bi);

bi_m1 = GetBinomials(m1);
mat1 = diag(1./bi_m1);

bi_m2 = GetBinomials(m2);
mat2 = diag(1./bi_m2);

fxy = mat1 * fxy_bi * mat2;

end