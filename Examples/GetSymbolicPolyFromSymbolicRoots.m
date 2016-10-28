function [f] = GetSymbolicPolyFromSymbolicRoots(root_mult_arr)
% Given the array of factors of f(x,y) 
%
% Inputs
%
% root_mult_arr : Factors of polynomial f(x,y)

syms x y

nDistinctFactors = size(root_mult_arr,1);

f = 1;

for i = 1:1:nDistinctFactors
    
    % Get the factor
    factor = root_mult_arr(i,1);
    
    % Get the multiplicity
    mult = root_mult_arr(i,2);
    
    temp_poly = factor.^mult;
    
    
    f = f * temp_poly; 
    
end



end