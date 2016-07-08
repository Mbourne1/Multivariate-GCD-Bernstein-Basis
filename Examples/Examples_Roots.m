function [fxy_matrix_exact,m] = Examples_Roots(ex_num)
% Examples_Roots(ex_num)
%
% Given an example number, return the matrix of coefficients of a
% polynomial f(x,y)
%
% Inputs.
%
% ex_num : Example Number



TYPE = 'From Coefficients';

switch TYPE 
    case 'From Roots'
        
        [fxy_matrix_exact,m] = Examples_Roots_FromRoots(ex_num);
        
    case 'From Coefficients'
        
        [fxy_matrix_exact,m] = Examples_Roots_FromCoefficients(ex_num);
        
end
    



end





