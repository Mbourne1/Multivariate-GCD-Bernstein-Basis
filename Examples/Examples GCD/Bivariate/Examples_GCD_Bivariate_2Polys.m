function [fxy,gxy,dxy,uxy,vxy,m,n,t] = Examples_GCD_Bivariate_2Polys(ex_num)
%
% % Inputs
%
% ex_num : Example Number
%
% % Outputs
%
% [fxy, gxy] : Coefficients of polynomial f(x,y) and g(x,y)
%
% dxy : Coefficients of polynomial d(x,y)
%
% [uxy, vxy] : Coefficients of polynomial u(x,y) and v(x,y)
%
% m : Degree of polynomial f(x,y)
%
% n : Degree of polynomial g(x,y)
%
% t : Degree of polynomial d(x,y)



EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
    case 'From Roots'
        [fxy, gxy, dxy, uxy, vxy, m, n, t] = Examples_GCD_FromRoots_Bivariate_2Polys(ex_num);
        
    case 'From Coefficients'
        [fxy, gxy, dxy, uxy, vxy, m, n, t] = Examples_GCD_FromCoefficients_Bivariate_2Polys(ex_num);
        
    otherwise
        error('err');
end


end

