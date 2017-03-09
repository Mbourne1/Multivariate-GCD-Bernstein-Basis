function [fxy, gxy, hxy, dxy, uxy, vxy, wxy, m, n, o, t] = Examples_GCD_FromCoefficients_Bivariate_3Polys(ex_num)
% Given an example number, get coefficient matrices of polynomials f(x,y),
% g(x,y) and their greatest common divisor d(x,y).
%
% Input. 
% 
% ex_num : (String) Example Number
%
% Output.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of polynomial h(x,y)
%
% dxy : (Matrix) Coefficients of polynomial d(x,y), the greatest common
% divisor of f(x,y), g(x,y) and h(x,y)
%
% uxy : (Matrix) Coefficients of polynomial u(x,y), the quotient polynomial
% such that f(x,y)/u(x,y) = d(x,y)
%
% vxy : (Matrix) Coefficients of polynomial u(x,y), the quotient polynomial
% such that g(x,y)/v(x,y) = d(x,y)
%
% wxy : (Matrix) Coefficients of polynomial w(x,y), the quotient polynomial
% such that h(x,y)/w(x,y) = d(x,y)

% Initialise symbolic variables
syms x y;

% Add Examples folder
addpath(genpath('../Examples'));

% Get the sets of factors and multiplicities of the polynomials f(x,y),
% g(x,y), h(x,y), d(x,y) u(x,y), v(x,y) and w(x,y)
[f_root_mult_arr, g_root_mult_arr, h_root_mult_arr, d_root_mult_arr,...
    u_root_mult_arr, v_root_mult_arr, w_root_mult_arr ] = GCD_Examples_Bivariate_3Polys(ex_num);

% Get coefficients of f(x,y), g(x,y) and h(x,y) from the symbolic factors
fxy = GetCoefficientsFromSymbolicRoots_Bivariate(f_root_mult_arr);
gxy = GetCoefficientsFromSymbolicRoots_Bivariate(g_root_mult_arr);
hxy = GetCoefficientsFromSymbolicRoots_Bivariate(h_root_mult_arr);

% Get coefficients of d(x,y) from symbolic factors
dxy = GetCoefficientsFromSymbolicRoots_Bivariate(d_root_mult_arr);

% Get coefficients of u(x,y), v(x,y) and w(x,y) from symbolic factors
uxy = GetCoefficientsFromSymbolicRoots_Bivariate(u_root_mult_arr);
vxy = GetCoefficientsFromSymbolicRoots_Bivariate(v_root_mult_arr);
wxy = GetCoefficientsFromSymbolicRoots_Bivariate(w_root_mult_arr);
             
% Get symbolic forms of the polynomials
symbolic_d = GetSymbolicPolyFromSymbolicRoots_Bivariate(d_root_mult_arr);

symbolic_f = GetSymbolicPolyFromSymbolicRoots_Bivariate(f_root_mult_arr);
symbolic_g = GetSymbolicPolyFromSymbolicRoots_Bivariate(g_root_mult_arr);
symbolic_h = GetSymbolicPolyFromSymbolicRoots_Bivariate(h_root_mult_arr);

% Get total degree of f(x,y), g(x,y), h(x,y) and d(x,y)
m = double(feval(symengine, 'degree', symbolic_f));
n = double(feval(symengine, 'degree', symbolic_g));
o = double(feval(symengine, 'degree', symbolic_h));
t = double(feval(symengine, 'degree', symbolic_d));

% Print coefficients to screen
% display(symbolic_f);
% display(symbolic_g);
% display(symbolic_h);
% display(symbolic_d);

end



