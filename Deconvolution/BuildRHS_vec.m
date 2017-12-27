
function [vRHS] = BuildRHS_vec(arr_fxy)
% Used in the batch deconvolution problem. This function constructs the
% right-hand side vector.
% 
% % Inputs
%
% arr_fxy : (Array of Matrices) Each matrix contains coefficients of a
% polynomial f_{i}(x,y)
%
% % Outputs
%
% vRHS : (Vector) Vector containing coefficients of the polynomials
% f_{0},..., f_{n-1} where n is the number of polynomials in the array.

% Get number of polynomials in the array 
nPolys_arr_fxy = length(arr_fxy);

% Initialise the RHS vector.
vRHS = [];

% For all but the final polynomial in the array, add the coefficients of
% f_{i}(x,y) to the RHS vector.
for i = 1 : 1 : nPolys_arr_fxy - 1
    
    fxy = arr_fxy{i};
    
    temp_vec = GetAsVector(fxy);
    
    vRHS = [vRHS ; temp_vec];
    
end





end