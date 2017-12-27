function fxy_matrix = GetAsMatrix_Version2(fxy_vec, m1, m2)
% Given the vector of coefficients of the polynomial f(x,y), 
% format the coefficients as a matrix. 
%
% % Inputs
%
% fxy_vec : (Vector) Coefficients of f(x,y) as a vector
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y
%
% % Outputs
%
% fxy_matrix : (Matrix) Coefficients of f(x,y) as a matrix.


% Initialise an empty matrix fxy
fxy_matrix = zeros(m1+1, m2+1);

for i2 = 0:1:m2
   
    start_position = (i2*(m1 + 1)) + 1;
    
    temp_vec = fxy_vec(start_position : start_position + m1);
    
    fxy_matrix(:,i2+1) = temp_vec;
    
end


end