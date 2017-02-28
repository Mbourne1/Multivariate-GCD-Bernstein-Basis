function [dxy] = GetGCD_Coefficients_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, t, t1, t2)
% 
%
% % Inputs
%
% [fxy, gxy, hxy] : Coefficients of the polynomial f(x,y), g(x,y) and
% h(x,y)
%
% [uxy, vxy, wxy] : Coefficients of the polynomial u(x,y), v(x,y) and
% w(x,y)
%
% [m, n, o] : Total degree of polynomial f(x,y), g(x,y) and h(x,y)
%
% t : Total degree of polynomial d(x,y) 
% 
% t1 : Degree of polynomial d(x,y) with respect to x
%
% t2 : Degree of polynomial d(x,y) with respect to y
%
% % Outputs
%
% dxy : Coefficients of the polynomial d(x,y)

global SETTINGS


switch SETTINGS.CALC_METHOD
    
    case 'Relative'
        
        dxy = GetGCD_Coefficients_Relative_3Polys(uxy, vxy, wxy, fxy, gxy, hxy, t1, t2);
        
    case 'Total'
        
        dxy = GetGCD_Coefficients_Total_3Polys(uxy, vxy, wxy, fxy, gxy, hxy, m, n, o, t);
        
    otherwise
        
        error('error')
        
        
end

end

