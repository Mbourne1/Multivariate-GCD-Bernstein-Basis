function [fxy,gxy,dxy,uxy,vxy,m,n,t] = Examples_GCD(ex_num)
%
% Inputs
%
% ex_num : Example Number
%
% Outputs
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial f(x,y)
%
% dxy : Coefficients of polynomial f(x,y)
%
% uxy : Coefficients of polynomial f(x,y)
%
% vxy : Coefficients of polynomial f(x,y)

EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
    case 'From Roots'
        [fxy,gxy,dxy,uxy,vxy,m,n,t] = Examples_GCD_FromRoots(ex_num);
        
    case 'From Coefficients'
        [fxy,gxy,dxy,uxy,vxy,m,n,t] = Examples_GCD_FromCoefficients(ex_num);
        
    otherwise
        error('err');
end


end

