function fxy_matrix = GetAsMatrix(fxy_vec, m1, m2)
% Given the vector of coefficients of the polynomial f(x,y),
% format the coefficients as a matrix.
%
% Note : This method inserts coefficients into a matrix in the order [0,0] [1,0] [0,1]
%
% % Inputs
%
% fxy_vec : (Vector) Coefficients of the polynomial f(x,y) written as a vector.
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y
%
% % Outputs
%
% fxy_matrix : (Matrix) Coefficients of the polynomial f(x,y) as a matrix.

global SETTINGS


switch SETTINGS.VECTORISATION_METHOD
    
    case 'Version 1'
        
        fxy_matrix = GetAsMatrix_Version1(fxy_vec, m1, m2);
        
    case 'Version 2'
        
        fxy_matrix = GetAsMatrix_Version2(fxy_vec, m1, m2);
        
    otherwise
        
        error('Error')
        
end
end