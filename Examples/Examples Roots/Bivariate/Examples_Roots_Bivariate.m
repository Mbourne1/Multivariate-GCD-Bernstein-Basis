function [fxy,m] = Examples_Roots_bivar(ex_num)
% Examples_Roots(ex_num)
%
% Given an example number, return the matrix of coefficients of a
% polynomial f(x,y)
%
% % Inputs.
%
% ex_num : Example Number as string
%
% % Outputs.
%
% fxy : Coefficients of matrix f(x,y)
%
% m : Total degree of f(x,y)




TYPE = 'From Coefficients';

switch TYPE 
    case 'From Roots'
        
        [fxy,m] = Examples_Roots_FromRoots(ex_num);
        
    case 'From Coefficients'
        
        [fxy,m] = Examples_Roots_FromCoefficients(ex_num);
        
end
    



end





