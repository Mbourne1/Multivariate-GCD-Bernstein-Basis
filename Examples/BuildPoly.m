function [fxy_bi] = BuildPoly(root_mult_matrix)
% Given a set of roots for f(x,y), build the polynomial so that the output
% is a vector of coefficients for the polynomial in the modified Bernstein 
% basis.

% Inputs
%
% root_mult_array : (Matrix) Consists of polynomial roots and corresponding
% multiplicity
% 
%
% Outputs.
%
% fxy_bi : (Matrix) Coefficients of polynomial f(x,y) in scaled Bernstein 
% form

% Calculate the number of distinct roots of the polynomial.
[nDistinctFactors, ~] = size(root_mult_matrix);

% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.

fxy_bi = 1;

for k = 1 : 1 : nDistinctFactors
    
    % Get the kth factors
    factor = root_mult_matrix(k,1);
    
    % Get the multiplicity of the kth factor
    mult = root_mult_matrix(k,2);
    
    % Convolve m times
    temp_poly = B_conv(factor, mult);
    
    % Build fxy_bi
    fxy_bi = conv(fxy_bi , temp_poly) ;
    
end


end