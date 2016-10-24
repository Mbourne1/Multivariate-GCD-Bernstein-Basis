function fx = GetCoefficients(arr_factors_fxy)
% Given an array of factors of f(x,y) get the coefficients of f(x,y) in the
% bivariate bernstein form.

syms x y;

% for each factor in the array of factors get coefficients in Bernstein form
for i = 1:1:length(arr_factors_fxy)
    
    % Get factor in power form
    power_factor = double(rot90(coeffs(arr_factors_fxy{i},[x,y],'All'),2));
    
    % Get the coefficients of factor in Bernstein form
    arr_factors{i} = PowerToBernstein_Bivariate(power_factor);
    
end

fx = arr_factors{1};

for i = 2:1:length(arr_factors_fxy)
    fx = Bernstein_Multiply_Bivariate(fx,arr_factors{i});
end


end