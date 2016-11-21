function [fx_Bern] = Power2Bern_univariate(fx_pwr)
% given a set of univariate coefficients in the power basis, obtain the
% coefficients in the Bernstein basis.

%       _
%    1 | |    
%    x | |
% x^{2}|_|

% Get the dimensions of the input polynomial
[r,c] = size(fx_pwr);

% if more columns than rows, then we have a row vector rather than a column
% vector
if (c > r) 
    fx_pwr = fx_pwr' ;
    prob_type = 'row_vector';
else 
    prob_type = 'col_vector';
end

% Get the dimesions of the input polynomial
[r,c] = size(fx_pwr);

% Get the degree of polynomial f
m = r - 1;

% Build the conversion Matrix

% for each row
for i = 0:1:m
    for j = 0:1:m
        if j>= i;
            mat(i+1,j+1) = nchoosek(m-i,j-i) ./ nchoosek(m,j);
        end
    end
end

% Alternative form
for i = 0:1:m
    for j = 0:1:m
        if i >= j;
            mat2(i+1,j+1) = nchoosek(i,j) ./ nchoosek(m,j);
        end
    end
end

% Obtain the Bernstein Coefficients by multiplying the power coefficients
% by the basis conversion matrix.
mat = mat'
fx_Bern = mat * fx_pwr ;

% If we started with a row vector
if strcmp(prob_type,'row_vector')
    % transpose so that output bernstein coefficients are also a row
    % vector.
    fx_Bern = fx_Bern';
end