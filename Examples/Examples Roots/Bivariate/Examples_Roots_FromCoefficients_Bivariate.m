function [fxy, m] = Examples_Roots_FromCoefficients_Bivariate(ex_num)
%
% % Inputs
% 
% ex_num : (String) Example number as a string
%
% % Outputs
%
% fxy : (Matrix) Coefficients of f(x,y)
%
% m : (Int) Total degree of f(x,y)

syms x y ;

% Add examples folder to current path
addpath(genpath('../Examples'));

% Get the factor and multiplicity array for polynomial f(x,y)
root_mult_arr = Roots_Examples_Bivariate(ex_num);

% return the coefficient matrix of f(x,y)
fxy = GetCoefficientsFromSymbolicRoots_Bivariate(root_mult_arr);

% Get the symbolic polynomial in power basis
symbolic_fxy = GetSymbolicPolyFromSymbolicRoots_Bivariate(root_mult_arr);

% Get the total degree of f(x,y)
m = double(feval(symengine, 'degree', symbolic_fxy));

fprintf('%s : Polynomial f(x,y). \n',mfilename);
display(symbolic_fxy);


end