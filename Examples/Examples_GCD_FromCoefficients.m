function [fxy,gxy,dxy,m,n,t] = Examples_GCD_FromCoefficients(ex_num)
% Given an example number, get coefficient matrices of polynomials f(x,y),
% g(x,y) and their greatest common divisor d(x,y).
%
% Input. 
% 
% ex_num

% Initialise symbolic variables
syms x y;

% Get example polynomials in symbolic form
switch ex_num
    case '1'
        
        d = (x^2 + y^2 + 0.95)^2 * (x-0.8624)^3;
        u = (x + 0.89547)^2;
        v = (y + 0.72531)^3;
        f = d*u;
        g = d*v;
        
    case '2'
        
        d = (x + y - 0.29574)^3 * (x-2.7246);
        u = (x+0.57432)^3 * (y-1.23145)^3;
        v = (y-0.547321)^2 * (x-0.7891011)^2;
        f = d*u;
        g = d*v;
        
    case '3'
        
        d = (x + 2*y +0.5678652)^3;
        u = (x-0.5678914)^2;
        v = (x-0.1234565)^2;
        f = u*d;
        g = v*d;
        
    case '4'
        
        d = (x + 2*y +0.5678652)^5;
        u = (x-0.5678914)^2;
        v = (x-0.1234565)^2;
        f = u*d;
        g = v*d;
        
        
    otherwise
        error('err')
        
end

% Display symbolic polynomials f and g
display(f)
display(g)

% Get degree of f(x,y), g(x,y) and d(x,y)
m = double(feval(symengine, 'degree', f));
n = double(feval(symengine, 'degree', g));
t = double(feval(symengine, 'degree', d));

% Get coefficients of f(x,y), g(x,y) and d(x,y)
fxy_pwr = double(rot90(coeffs(f,[x,y],'All'),2));
gxy_pwr = double(rot90(coeffs(g,[x,y],'All'),2));
dxy_pwr = double(rot90(coeffs(d,[x,y],'All'),2));

% Convert coefficients to Bernstein form.
fxy = Power2Bernstein_Bivariate(fxy_pwr);
gxy = Power2Bernstein_Bivariate(gxy_pwr);
dxy = Power2Bernstein_Bivariate(dxy_pwr);

end
