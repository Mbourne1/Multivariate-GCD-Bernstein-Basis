function [Bern_coeff] = Power2Bern_Bivariate(fxy_pwr)
% given a set of univariate coefficients in the power basis, obtain the
% coefficients in the Bernstein basis.

%       1 y y^2
%       _ _ _
%    1 |_|_|_| 
%    x |_|_|_|
% x^{2}|_|_|_|

[r,c] = size(fxy_pwr);

% Get the degree of polynomial f
m = r - 1;
n = c - 1;

% Build the conversion Matrix

% for each row
for i = 0:1:m
    for j = 0:1:m
        if j>= i;
            mat1(i+1,j+1) = nchoosek(m-i,j-i) ./ nchoosek(m,j);
        end
    end
end


% for each col
for i = 0:1:n
    for j = 0:1:n
        if j>= i;
            mat2(i+1,j+1) = nchoosek(n-i,j-i) ./ nchoosek(n,j);
        end
    end
end


Bern_coeff = mat1' * fxy_pwr * mat2


