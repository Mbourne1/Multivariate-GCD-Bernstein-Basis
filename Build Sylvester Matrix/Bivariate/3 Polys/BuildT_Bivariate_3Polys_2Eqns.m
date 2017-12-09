function T_k1k2 = BuildT_Bivariate_3Polys_2Eqns(fxy, gxy, hxy, k1, k2)
%
% % Inputs
%
% fxy : (Matrix)
%
% gxy : (Matrix)
%
% hxy : (Matrix)
%
% k1 : (Int)
%
% k2 : (Int)

% Get degree of polynomial f(x,y) and g(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);


T1 = BuildT1_Bivariate(fxy, n1 - k1, n2 - k2);
T2 = BuildT1_Bivariate(gxy, m1 - k1, m2 - k2);
T3 = BuildT1_Bivariate(fxy, o1 - k1, o2 - k2);
T4 = BuildT1_Bivariate(hxy, m1 - k1, m2 - k2);

T_k1k2 = [blkdiag(T1,T3) [T2; T4]];

end