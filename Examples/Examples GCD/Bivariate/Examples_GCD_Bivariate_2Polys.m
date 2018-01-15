function [fxy,gxy,dxy,uxy,vxy,m,n,t] = Examples_GCD_Bivariate_2Polys(ex_num)
%
% % Inputs
%
% ex_num : Example Number
%
% % Outputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% dxy : (Matrix) Coefficients of polynomial d(x,y)
%
% uxy : (Matrix) Coefficients of polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of polynomial v(x,y)
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)
%
% t : (Int) Total degree of polynomial d(x,y)




[fxy, gxy, dxy, uxy, vxy, m, n, t] = ...
    Examples_GCD_FromCoefficients_Bivariate_2Polys(ex_num);



end

