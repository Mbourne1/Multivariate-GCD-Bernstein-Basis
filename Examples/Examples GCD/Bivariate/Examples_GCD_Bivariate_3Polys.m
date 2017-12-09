function [fxy, gxy, hxy, dxy, uxy, vxy, wxy, m, n, o, t] = Examples_GCD_Bivariate_3Polys(ex_num)
%
% Inputs
%
% ex_num : (String) Example Number
%
% Outputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of polynomial h(x,y)
%
% dxy : (Matrix) Coefficients of polynomial d(x,y), the greatest common
% divisor of f(x,y) g(x,y) and h(x,y)
%
% uxy : (Matrix) Coefficients of polynomial u(x,y), the quotient
% polynomial such that f(x,y)/u(x,y) = d(x,y)
%
% vxy : (Matrix) Coefficients of polynomial v(x,y), the quotient polynomial
% such that g(x,y)/v(x,y) = d(x,y)
%
% wxy : (Matrix) Coefficients of polynomial w(x,y), the quotient polynomial
% such that h(x,y)/w(x,y) = d(x,y)


[fxy, gxy, hxy, dxy, uxy, vxy, wxy, m, n, o, t] = Examples_GCD_FromCoefficients_Bivariate_3Polys(ex_num);



end

