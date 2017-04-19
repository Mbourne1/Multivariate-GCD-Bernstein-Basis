function [DTQ] = BuildDTQ_Bivariate_2Polys_Version2(fxy, gxy, k1, k2)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% k1 : (Int)
%
% k2 : (Int)


% Get degree of f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of g(x,y) with respect to x and y
[n1, n2] = GetDegree_Bivariate(gxy);

% Build matrices D^{-1}T(f)Q and D^{-1}T(g)Q
DT1Q1 = BuildDT1Q1_Bivariate_Version2(fxy, n1-k1, n2-k2);
DT2Q2 = BuildDT1Q1_Bivariate_Version2(gxy, m1-k1, m2-k2);

% Build DTQ
DTQ = [DT1Q1 DT2Q2];



end

