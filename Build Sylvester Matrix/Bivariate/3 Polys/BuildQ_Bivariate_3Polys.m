function Q = BuildQ_Bivariate_3Polys(m1, m2, n1, n2, o1, o2, k1, k2)
% Build the block diagonal matrix \tilde{Q}_{k_{1}, k_{2}}
%
% % Inputs
%
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


Q_ok = BuildQ1_Bivariate(o1 - k1, o2 - k2);

Q_nk = BuildQ1_Bivariate(n1 - k1, n2 - k2);

Q_mk = BuildQ1_Bivariate(m1 - k1, m2 - k2);

Q = blkdiag(Q_nk, Q_ok, Q_mk);


end