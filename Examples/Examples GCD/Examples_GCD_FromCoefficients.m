function [fxy,gxy,dxy,uxy,vxy,m,n,t] = Examples_GCD_FromCoefficients(ex_num)
% Given an example number, get coefficient matrices of polynomials f(x,y),
% g(x,y) and their greatest common divisor d(x,y).
%
% Input. 
% 
% ex_num : Example Number
%
% Output.
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% dxy : Coefficients of polynomial d(x,y)
%
% uxy : Coefficients of polynomial u(x,y)
%
% vxy : Coefficients of polynomial v(x,y)



% Initialise symbolic variables
syms x y;

addpath('../Examples')
[f_root_mult_arr,g_root_mult_arr,d_root_mult_arr,...
    u_root_mult_arr,v_root_mult_arr] = Bivariate_GCD_Examples(ex_num);


fxy = GetCoefficientsFromSymbolicRoots(f_root_mult_arr);
gxy = GetCoefficientsFromSymbolicRoots(g_root_mult_arr);
dxy = GetCoefficientsFromSymbolicRoots(d_root_mult_arr);
uxy = GetCoefficientsFromSymbolicRoots(u_root_mult_arr);
vxy = GetCoefficientsFromSymbolicRoots(v_root_mult_arr);
      
             
symbolic_d = GetSymbolicPolyFromSymbolicRoots(d_root_mult_arr);
symbolic_f = GetSymbolicPolyFromSymbolicRoots(f_root_mult_arr);
symbolic_g = GetSymbolicPolyFromSymbolicRoots(g_root_mult_arr);

% Get degree of f(x,y), g(x,y) and d(x,y)
m = double(feval(symengine, 'degree', symbolic_f));
n = double(feval(symengine, 'degree', symbolic_g));
t = double(feval(symengine, 'degree', symbolic_d));

display(symbolic_f);
display(symbolic_g);
display(symbolic_d);

end



