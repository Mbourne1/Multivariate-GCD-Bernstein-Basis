function [fxy_matrix] = GetAsMatrix_Version1(fxy_vec, m1, m2)
%
% % Inputs
%
% fxy_vec : (Vector) Coefficients of the polynomial f(x,y)
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y
%
% % Outputs
%
% fxy_matrix : (Matrix) Coefficients of polynomial f(x,y)
 

% Initialise an empty matrix fxy
fxy_matrix = zeros(m1+1, m2+1);

% Intialise a counter which will go through each entry of f_vec (The vector
% of coefficients of of f).
count = 1;

% Get number of diagonals in the matrix fxy.
nDiagonals = (m1+1) + (m2+1) -1;


for tot = 0 : 1 : nDiagonals - 1
    for i = tot : -1 : 0

        j = tot-i;
        
        if i > m1 || j> m2
            
            % restrict to only the i and j values within the matrix.
            
        else
            
            fxy_matrix(i + 1, j + 1) = fxy_vec(count);
            count = count + 1;
        end
    end
end


