function D = BuildD_Bivariate_3Polys_3Eqns(m1, m2, n1, n2, o1, o2, k1, k2)
%
% % Inputs
%
% m1 : (Int)
%
% m2 : (Int)
%
% n1 : (Int)
%
% n2 : (Int)
%
% o1 : (Int)
%
% o2 : (Int)
%
% k1 : (Int)
%
% k2 : (Int)


D1 = BuildD_Bivariate_2Polys(m1, m2, n1 - k1, n2 - k2);
D2 = BuildD_Bivariate_2Polys(m1, m2, o1 - k1, o2 - k2);
D3 = BuildD_Bivariate_2Polys(n1, n2, o1 - k1, o2 - k2);

D = blkdiag(D1, D2, D3);

end