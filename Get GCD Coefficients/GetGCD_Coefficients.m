function [dxy] = GetGCD_Coefficients(fxy,gxy,uxy,vxy,m,n,t,t1,t2)
% 
%
% % Inputs
%
% fxy : Coefficients of the polynomial f(x,y)
%
% gxy : Coefficients of the polynomial g(x,y)
%
% uxy : Coefficients of the polynomial u(x,y)
%
% vxy : Coefficients of the polynomial v(x,y)
%
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomial g(x,y)
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
        
        dxy = GetGCD_Coefficients_Relative(uxy,vxy,fxy,gxy,t1,t2);
        
    case 'Total'
        
        dxy = GetGCD_Coefficients_Total(uxy,vxy,fxy,gxy,m,n,t);
        
    otherwise
        
        error('error')
        
        
end

end

