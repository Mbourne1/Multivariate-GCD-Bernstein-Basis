function [fxy,m] = Examples_Roots_FromCoefficients(ex_num)


syms x y ;

switch ex_num
    case '1'
        f =  (x-0.5)^2 * (x-0.789) * (y-0.5);
    otherwise
        error([mfilename ' : ' sprintf('%s is not a valid example number',ex_num)])
end
            
display f

% Get total degree of f
m = double(feval(symengine, 'degree', f));

% Get coefficients of f(x,y) in power basis
fxy_pwr = double(rot90(coeffs(f,[x,y],'All'),2));

% Get coefficients of f(x,y) in Bernstein basis
fxy = Power2Bernstein_Bivariate(fxy_pwr);

end