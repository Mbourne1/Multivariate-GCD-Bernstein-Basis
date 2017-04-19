function fxy = GetCoefficients_Bivariate(arr_factors_fxy)
% Given an array of factors of f(x,y) get the coefficients of f(x,y) in the
% bivariate bernstein form.
%
% % Inputs
% 
% arr_factors_fxy : Flat array of factors of f(x,y). Factors of 
% multiplicity 'm' are repeated 'm' times in this array.
%
%
% % Outputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)


syms x y;

% Get number of factors 
nFactors = size(arr_factors_fxy,1);

% For each factor in the array of factors get coefficients in Bernstein form
arr_factors = cell(nFactors,1);
for i = 1:1:nFactors
    
    % Get factor in power form
    power_factor = double(rot90(coeffs(arr_factors_fxy{i},[x,y],'All'),2));
    
    % Get the coefficients of factor in Bernstein form
    arr_factors{i} = PowerToBernstein_Bivariate(power_factor);
    
end


temp_coefficient_matrix = arr_factors{1};


for i = 2:1:nFactors
    temp_coefficient_matrix = Bernstein_Multiply_Bivariate(temp_coefficient_matrix, arr_factors{i});
end

fxy = temp_coefficient_matrix;

end