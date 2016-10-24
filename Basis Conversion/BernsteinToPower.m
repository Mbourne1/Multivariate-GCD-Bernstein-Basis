function [fx_Pwr] = BernsteinToPower(fx_Brn)
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

% Get the degree of input polynomial f(x) in Bernstein basis.
m = size(fx_Brn,1) - 1;


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


end