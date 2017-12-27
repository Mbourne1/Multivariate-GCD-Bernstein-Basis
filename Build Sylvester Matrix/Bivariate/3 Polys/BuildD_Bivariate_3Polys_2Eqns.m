function D = BuildD_Bivariate_3Polys_2Eqns(m1, m2, n1, n2, o1, o2, k1, k2)
%
% % Inputs
%
% m1 : (Int) Degree of f(x, y) with respect to x
%
% m2 : (Int) Degree of f(x, y) with respect to y
%
% n1 : (Int) Degree of g(x, y) with respect to x
%
% n2 : (Int) Degree of g(x, y) with respect to y
%
% o1 : (Int) Degree of h(x, y) with respect to x
%
% o2 : (Int) Degree of h(x, y) with respect to y
%
% k1 : (Int) Index of subresultant matrix
%
% k2 : (Int) Index of subresultant matrix


D1 = BuildD_Bivariate_2Polys(m1, m2, n1 - k1, n2 - k2);

D2 = BuildD_Bivariate_2Polys(m1, m2, o1 - k1, o2 - k2);

% Build block diagonal matrix of D1 and D2
D = blkdiag(D1,D2);

end