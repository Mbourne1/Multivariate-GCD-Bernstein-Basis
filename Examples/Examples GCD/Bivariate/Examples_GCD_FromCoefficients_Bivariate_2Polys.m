function [fxy, gxy, dxy, uxy, vxy, m, n, t] = Examples_GCD_FromCoefficients_Bivariate_2Polys(ex_num)
% Given an example number, get coefficient matrices of polynomials f(x,y),
% g(x,y) and their greatest common divisor d(x,y).
%
% Input.
%
% ex_num : Example Number
%
% Output.
%
% [fxy, gxy] : Coefficients of polynomial f(x,y) and g(x,y)
%
% dxy : Coefficients of polynomial d(x,y)
%
% [uxy, vxy] : Coefficients of polynomial u(x,y) and v(x,y)



% Initialise symbolic variables
syms x y;

addpath(genpath('../Examples'));

[f_root_mult_arr,g_root_mult_arr,d_root_mult_arr,...
    u_root_mult_arr,v_root_mult_arr] = GCD_Examples_Bivariate_2Polys(ex_num);


% Get coefficients of polynomials f(x,y) and g(x,y)
fxy = GetCoefficientsFromSymbolicRoots_Bivariate(f_root_mult_arr);
gxy = GetCoefficientsFromSymbolicRoots_Bivariate(g_root_mult_arr);

% Get coefficients of polynomial d(x,y)
dxy = GetCoefficientsFromSymbolicRoots_Bivariate(d_root_mult_arr);

% Get coefficients of the polynomials u(x,y) and v(x,y)
uxy = GetCoefficientsFromSymbolicRoots_Bivariate(u_root_mult_arr);
vxy = GetCoefficientsFromSymbolicRoots_Bivariate(v_root_mult_arr);


% Get Symbolic polynomials f(x,y), g(x,y) and d(x,y)
symbolic_d = GetSymbolicPolyFromSymbolicRoots_Bivariate(d_root_mult_arr);
symbolic_f = GetSymbolicPolyFromSymbolicRoots_Bivariate(f_root_mult_arr);
symbolic_g = GetSymbolicPolyFromSymbolicRoots_Bivariate(g_root_mult_arr);

% Get degree of f(x,y), g(x,y) and d(x,y)
m = double(feval(symengine, 'degree', symbolic_f));
n = double(feval(symengine, 'degree', symbolic_g));
t = double(feval(symengine, 'degree', symbolic_d));

display(symbolic_f);
display(symbolic_g);
display(symbolic_d);

end



