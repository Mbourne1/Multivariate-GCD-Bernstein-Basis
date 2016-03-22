function [f_bi] = BuildPoly(root_mult_matrix)
% Given a set of roots for f(x,y), build the polynomial so that the output
% is a vector of coefficients for the polynomial in the modified Bernstein 
% basis.

% Calculate the number of distinct roots of the polynomial.
[nEntries,~] = size(root_mult_matrix);

% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.

f_bi = 1;
for k = 1:1:nEntries
    w = B_conv(root_mult_matrix(k,1),root_mult_matrix(k,2));
    f_bi = conv(f_bi,w) ;
end


end