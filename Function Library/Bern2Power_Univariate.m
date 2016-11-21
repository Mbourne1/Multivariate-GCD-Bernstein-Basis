function [fx_Pwr] = Bern2Power_Univariate(fx_Brn)
% Given a column vector of coefficients of the polynomial f(x) in Bernstein
% form. Get the column vector of coefficients of the same polynomial in 
% power form.
%
% % Input:
%
% fx_Brn : Column vector of coefficients of polynomial f(x) in Bernstein form.
%                __
%     B^{m}_{0} |__|    
%     B^{m}_{1} |__|
%       ...     |__|
%     B^{m}_{m} |__|
%
% % Outputs : 
%
% fx_Pwr : Column vector of coefficients of polynomial f(x) in Power form.
%
%

% Get the dimensions of the input polynomial
[r,c] = size(fx_Brn);

% if more columns than rows, then we have a row vector rather than a column
% vector
if (c > r) 
    fx_Brn = fx_Brn' ;
    prob_type = 'row_vector';
else 
    prob_type = 'col_vector';
end

% Get the dimesions of the input polynomial
[r,~] = size(fx_Brn);

% Get the degree of polynomial f
m = r - 1;

% Build the conversion Matrix
mat = zeros(m+1,m+1);
% for each row
for i = 0:1:m
    for j = 0:1:m
        if j>= i;
            mat(i+1,j+1) = nchoosek(m-i,j-i) ./ nchoosek(m,j);
        end
    end
end

% Obtain the Bernstein Coefficients by multiplying the power coefficients
% by the basis conversion matrix.
mat = pinv(mat');
fx_Pwr =  mat * fx_Brn  ;

% If we started with a row vector
if strcmp(prob_type,'row_vector')
    % transpose so that output bernstein coefficients are also a row
    % vector.
    fx_Pwr = fx_Pwr';
end