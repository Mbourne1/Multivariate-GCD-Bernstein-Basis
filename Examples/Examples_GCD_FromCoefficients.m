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
        
        
        f = (x^2 + y^2 + 1)^2 * (x+1) * (y-1) * (y-3);
        g = (x^2 + y^2 + 1)^2 * (x+1) * (y-2);
        d = (x^2 + y^2 + 1)^2 * (x+1);
        
    case '2'
        
        f = (x+y+1) * (x+2) * (y-1);
        g = (x+y+1) * (x+1) * (y-2);
        d = (x+y+1);
    case '3'
        f = (x + 1) * (x+y+1) * (y-0.2);
        g = (x + 1) * (x+y+1) * (y-0.4) * (x-0.3);
        d = (x + 1) * (x+y+1);
    case '4'
        f = (x + 1) * (x+y+1)^2 * (y-0.2);
        g = (x + 1) * (x+y+1)^2 * (y-0.4) * (x-0.3);
        d = (x + 1) * (x+y+1)^2;
    otherwise
        error('err')
        
end

% Display symbolic polynomials f and g
display f
display g

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
