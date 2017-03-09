function [fxy] = GetCoefficientsFromSymbolicRoots_Bivariate(root_mult_arr_f)
% Given the set of symbolic factors and corresponding multiplicities, get
% the matrix containing coefficients of the polynomial f(x,y)
%
% % Inputs
%
% root_mult_arr_f : (Matrix) of symbolic factors f_{i} and corresponding
% multiplicities m_{i}.
%
% % Outputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)


% Get array of symbolic factors of the polynomial
sym_factor_arr_f = GetFactors_Bivariate(root_mult_arr_f);

% Get the coefficients of f(x,y) in Bernstein form
fxy = GetCoefficients_Bivariate(sym_factor_arr_f);

end