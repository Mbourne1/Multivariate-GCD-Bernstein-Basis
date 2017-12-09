function T_k1k2 = BuildT_Bivariate_3Polys_3Eqns(fxy, gxy, hxy, k1, k2)
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


%T = [T1 T2 T3 ; T4 T5 T6; T7 T8 T9]


T1 = BuildT1_Bivariate(fxy, n1 - k1, n2 - k2);
T3 = BuildT1_Bivariate(gxy, m1 - k1, m2 - k2);
T5 = BuildT1_Bivariate(fxy, o1 - k1, o2 - k2);
T6 = BuildT1_Bivariate(hxy, m1 - k1, m2 - k2);
T7 = BuildT1_Bivariate(hxy, n1 - k1, n2 - k2);
T8 = BuildT1_Bivariate(gxy, o1 - k1, o2 - k2);

% Construct the zero partitions
nRows = (m1 + n1 - k1 + 1) * (m2 + n2 - k2 + 1);
nCols = (o1 - k1 + 1) * (o2 - k2 + 1);
T2 = zeros(nRows, nCols);

nRows = (m1 + o1 - k1 + 1) * (m2 + o2 - k2 + 1);
nCols = (n1 - k1 + 1) * (n2 - k2 + 1);
T4 = zeros(nRows, nCols);


nRows = (n1 + o1 - k1 + 1) * (n2 + o2 - k2 + 1);
nCols = (m1 - k1 + 1) * (m2 - k2 + 1);
T9 = zeros(nRows, nCols);

T_k1k2 = [T1 T2 T3 ; T4 T5 T6; T7 -T8 T9];


end