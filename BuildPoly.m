function [f_bi] = BuildPoly(A)

% Calculate the number of distinct roots of the polynomial.
[r,c] = size(A);

% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.

f_bi = 1;
for k = 1:1:r
    w = B_conv(A(k,1),A(k,2));
    f_bi = conv(f_bi,w) ;
end


end