function fxy = GetCoefficients_Bivariate(arr_factors_fxy)
% Given an array of factors of f(x,y) get the coefficients of f(x,y) in the
% bivariate bernstein form.

syms x y;

nFactors = size(arr_factors_fxy,1);

% for each factor in the array of factors get coefficients in Bernstein form
for i = 1:1:nFactors
    
    % Get factor in power form
    power_factor = double(rot90(coeffs(arr_factors_fxy{i},[x,y],'All'),2));
    
    % Get the coefficients of factor in Bernstein form
    arr_factors{i} = PowerToBernstein_Bivariate(power_factor);
    
end

temp_coef_mat = arr_factors{1};

for i = 2:1:nFactors
    temp_coef_mat = Bernstein_Multiply_Bivariate(temp_coef_mat,arr_factors{i});
end

fxy = temp_coef_mat;

end