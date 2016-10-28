
function [arr_factors] = GetFactors(root_mult_arr)

% Given the set of factors and multiplicities, get an array of the factors
% where each factor is repeated m times.

arr_factors = {};

% Get number of unique factors
nUniqueFactors = size(root_mult_arr,1);



% for each row
for i = 1:1:nUniqueFactors
    
    % Get the symbolic factor
    sym_factor = root_mult_arr(i,1);
    
    % Get the multiplicity of the factor
    mult = root_mult_arr(i,2);
    
    arr_factors = [arr_factors ; repmat({sym_factor},1,mult)'];
    
end



end