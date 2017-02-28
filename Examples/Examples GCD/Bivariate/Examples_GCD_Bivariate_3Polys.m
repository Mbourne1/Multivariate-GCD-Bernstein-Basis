function [fxy, gxy, hxy, dxy, uxy, vxy, wxy, m, n, o, t] = Examples_GCD_Bivariate_3Polys(ex_num)
%
% Inputs
%
% ex_num : Example Number
%
% Outputs
%
% [fxy, gxy, hxy] : Coefficients of polynomial f(x,y), g(x,y) and h(x,y)
%
% dxy : Coefficients of polynomial f(x,y)
%
% [uxy, vxy, wxy] : Coefficients of polynomial u(x,y), v(x,y) and w(x,y)

EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
    case 'From Roots'
        
        [fxy, gxy, hxy, dxy, uxy, vxy, wxy, m, n, o, t] = Examples_GCD_FromRoots_3Polys(ex_num);
        
    case 'From Coefficients'
        
        [fxy, gxy, hxy, dxy, uxy, vxy, wxy, m, n, o, t] = Examples_GCD_FromCoefficients_Bivariate_3Polys(ex_num);
        
    otherwise
        error('err');
end


end

