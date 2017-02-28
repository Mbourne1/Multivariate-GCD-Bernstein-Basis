function [fxy_bi] = BuildPoly(root_mult_matrix)
% Given a set of roots for f(x,y), build the polynomial so that the output
% is a vector of coefficients for the polynomial in the modified Bernstein 
% basis.

% Inputs
%
% root_mult_array
%
% root : factor 
%
% Outputs.
%
% fxy_bi : Matrix of coefficients of polynomial f(x,y) in scaled Bernstein 
% form, ie. in the form 'a_{i,j} \binom{m1}{i}\binom{m2}{j}'

% Calculate the number of distinct roots of the polynomial.
[nEntries,~] = size(root_mult_matrix);

% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.

fxy_bi = 1;
for k = 1:1:nEntries
    
    root = root_mult_matrix(k,1);
    mult = root_mult_matrix(k,2);
    
    temp_poly = B_conv(root,mult);
    
    fxy_bi = conv(fxy_bi , temp_poly) ;
end


end