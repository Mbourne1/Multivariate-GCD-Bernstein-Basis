function [f] = GetSymbolicPolyFromSymbolicRoots_Bivariate(root_mult_arr)
% Given the array of symbolic factors of f(x,y) and corresponding 
% multiplicities return a symbolic polynomial 
%
% % Inputs
%
% root_mult_arr : (Matrix) Factors and multiplicities of polynomial f(x,y)
%
% % Outputs
%
% f : (Symbolic) Symbolic polynomial f(x,y)

syms x y

% Get number of distinct factors in f(x,y)
nDistinctFactors = size(root_mult_arr,1);

% Initialise symbolic polynomial f
f = 1;

for i = 1:1:nDistinctFactors
    
    % Get the ith factor
    factor = root_mult_arr(i,1);
    
    % Get the multiplicity of ith factor
    mult = root_mult_arr(i,2);
    
    % Build the symbolic polynomial f(x,y)
    temp_poly = factor.^mult;
    f = f * temp_poly; 
    
end



end