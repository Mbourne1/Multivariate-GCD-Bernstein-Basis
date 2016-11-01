function [fxy,m] = Examples_Roots_FromCoefficients(ex_num)


syms x y ;

addpath '../Examples'

% Get the factor and multiplicity array for polynomial f(x,y)
root_mult_arr = Bivariate_Roots_Examples(ex_num);

% return the coefficient matrix of f(x,y)
fxy = GetCoefficientsFromSymbolicRoots(root_mult_arr);

% Get the symbolic polynomial in power basis
symbolic_fxy = GetSymbolicPolyFromSymbolicRoots(root_mult_arr);

% Get the total degree of f(x,y)
m = double(feval(symengine, 'degree', symbolic_fxy));

fprintf('%s : Polynomial f(x,y). \n',mfilename);
display(symbolic_fxy);


end