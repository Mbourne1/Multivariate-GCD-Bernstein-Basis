function [arr_factors] = GetFactors_Bivariate(root_mult_arr)
%
% Given the set of factors and multiplicities, get an array of the factors
% where each factor is repeated m times.
%
% % Inputs
%
% root_mult_arr : (Matrix) containing symbolic factors of f(x,y) and their
% corresponding multiplicity
%
% % Outputs
%
% Array of factors f_{i} each repeated m_{i} times

% Initialise an array to store the factors
arr_factors = {};

% Get number of unique factors
nUniqueFactors = size(root_mult_arr,1);



% For each unique factor
for i = 1:1:nUniqueFactors
    
    % Get the symbolic factor
    sym_factor = root_mult_arr(i,1);
    
    % Get the multiplicity of the factor
    mult = root_mult_arr(i,2);
    
    % Insert the symbolic factor into the array m times
    arr_factors = [arr_factors ; repmat({sym_factor},1,mult)'];
    
end



end