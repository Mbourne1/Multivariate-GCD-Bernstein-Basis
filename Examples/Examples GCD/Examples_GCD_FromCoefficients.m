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
[f,g,d,u,v] = Bivariate_GCD_Examples(ex_num);


fxy = GetCoefficients(f);
gxy = GetCoefficients(g);
dxy = GetCoefficients(d);
uxy = GetCoefficients(u);
vxy = GetCoefficients(v);


symbolic_d = GetSymbolicPoly(d);
symbolic_f = GetSymbolicPoly(f);
symbolic_g = GetSymbolicPoly(g);

% Get degree of f(x,y), g(x,y) and d(x,y)
m = double(feval(symengine, 'degree', symbolic_f));
n = double(feval(symengine, 'degree', symbolic_g));
t = double(feval(symengine, 'degree', symbolic_d));

display(symbolic_f);
display(symbolic_g);
display(symbolic_d);

end



