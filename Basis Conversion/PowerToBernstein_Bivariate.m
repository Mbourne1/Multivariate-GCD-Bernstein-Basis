function [fxy_Bb] = PowerToBernstein_Bivariate(fxy_matrix)
% Given the coefficients of a bivariate polynomial f(x,y) in power basis,
% return the coefficients of the same polynomial in Bernstein basis.
%
% > Power2Bernstein_Bivariate([1 1; 1 0])

% Get the degree of polynomial f(x,y) with respect to x (m1) and with
% respect to y (m2)
[m1, m2] = GetDegree_Bivariate(fxy_matrix);


% Build the matrices A and B such that A * f(x,y) * B gives the
% coefficients of the polynomial in the Bernstein Basis.
A = zeros(m1+1, m1+1);
B = zeros(m2+1, m2+1);

for i = 0:1:m1
    for j = 0:1:m1
        if i >= j
            A(i+1,j+1) = nchoosek(i,j) ./ nchoosek(m1,j);
        end
    end
end

for i = 0:1:m2
    for j = 0:1:m2
        if i >=j
            B(i+1,j+1) = nchoosek(i,j) ./ nchoosek(m2,j);
        end
    end
end

fxy_Bb = A* fxy_matrix * B';



end