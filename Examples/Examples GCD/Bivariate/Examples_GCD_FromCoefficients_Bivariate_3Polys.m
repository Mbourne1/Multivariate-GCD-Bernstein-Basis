function [fxy, gxy, hxy, dxy, uxy, vxy, wxy, m, n, o, t] = Examples_GCD_FromCoefficients_Bivariate_3Polys(ex_num)
% Given an example number, get coefficient matrices of polynomials f(x,y),
% g(x,y) and their greatest common divisor d(x,y).
%
% Input. 
% 
% ex_num : Example Number
%
% Output.
%
% [fxy, gxy, hxy] : Coefficients of polynomial f(x,y), g(x,y) and h(x,y)
%
% dxy : Coefficients of polynomial d(x,y)
%
% [uxy, vxy, wxy] : Coefficients of polynomial u(x,y), v(x,y) and w(x,y)

% Initialise symbolic variables
syms x y;

addpath(genpath('../Examples'));

[f_root_mult_arr, g_root_mult_arr, h_root_mult_arr, d_root_mult_arr,...
    u_root_mult_arr, v_root_mult_arr, w_root_mult_arr ] = GCD_Examples_Bivariate_3Polys(ex_num);


fxy = GetCoefficientsFromSymbolicRoots_Bivariate(f_root_mult_arr);
gxy = GetCoefficientsFromSymbolicRoots_Bivariate(g_root_mult_arr);
hxy = GetCoefficientsFromSymbolicRoots_Bivariate(h_root_mult_arr);

dxy = GetCoefficientsFromSymbolicRoots_Bivariate(d_root_mult_arr);


uxy = GetCoefficientsFromSymbolicRoots_Bivariate(u_root_mult_arr);
vxy = GetCoefficientsFromSymbolicRoots_Bivariate(v_root_mult_arr);
wxy = GetCoefficientsFromSymbolicRoots_Bivariate(w_root_mult_arr);
             
symbolic_d = GetSymbolicPolyFromSymbolicRoots_Bivariate(d_root_mult_arr);

symbolic_f = GetSymbolicPolyFromSymbolicRoots_Bivariate(f_root_mult_arr);
symbolic_g = GetSymbolicPolyFromSymbolicRoots_Bivariate(g_root_mult_arr);
symbolic_h = GetSymbolicPolyFromSymbolicRoots_Bivariate(h_root_mult_arr);

% Get degree of f(x,y), g(x,y) and d(x,y)
m = double(feval(symengine, 'degree', symbolic_f));
n = double(feval(symengine, 'degree', symbolic_g));
o = double(feval(symengine, 'degree', symbolic_h));
t = double(feval(symengine, 'degree', symbolic_d));


display(symbolic_f);
display(symbolic_g);
display(symbolic_h);
display(symbolic_d);

end



