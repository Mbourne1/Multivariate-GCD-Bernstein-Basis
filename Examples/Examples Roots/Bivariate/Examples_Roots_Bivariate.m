function [fxy, m] = Examples_Roots_Bivariate(ex_num)
% Examples_Roots_Bivariate(ex_num)
%
% Given an example number, return the matrix of coefficients of a
% polynomial f(x,y)
%
% % Inputs.
%
% ex_num : (String) Example Number as string
%
% % Outputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% m : (Int) Total degree of f(x,y)


TYPE = 'From Coefficients';

switch TYPE 
    case 'From Roots'
        
        [fxy,m] = Examples_Roots_FromRoots_Bivariate(ex_num);
        
    case 'From Coefficients'
        
        [fxy,m] = Examples_Roots_FromCoefficients_Bivariate(ex_num);
        
end
    



end





