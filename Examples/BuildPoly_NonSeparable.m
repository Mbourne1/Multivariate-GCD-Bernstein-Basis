function [fxy] = BuildPoly_NonSeparable(factors_fxy_array)
% Given a set of roots for f build the polynomial so that we obtain a set
% of coefficients for the polynomial in the Bernstein basis.
% 
% Inputs.
%
% factors_fxy_array : (Array) Array of factors of f(x,y) 
%
% Outputs.
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)


% Get the number of distinct roots.
[nFactors,~] = size(factors_fxy_array);

% Initialise the product
temp_prod = factors_fxy_array{1,1};

% For every other root 2,...,
for i = 2:1:nFactors
    
   % Get coefficients of the factor as a vector
   temp_factor = factors_fxy_array{i,1};
   
   % Multiply product by the new factor
   temp_prod = Bernstein_Multiply_Bivariate(temp_prod, temp_factor ) ;
   
end

% Output the product of all roots : f(x,y) 
fxy = temp_prod;


end