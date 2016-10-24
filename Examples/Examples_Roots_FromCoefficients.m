function [fxy,m] = Examples_Roots_FromCoefficients(ex_num)


syms x y ;

switch ex_num
    case '1'
        sym_factor_mult_arr_f = [
            (x-0.5)     2 
            (x-0.789)   1
            (y-0.5)     1
            ];
        
        sym_factor_arr_f = GetFactors(sym_factor_mult_arr_f);
        
    otherwise
        error([mfilename ' : ' sprintf('%s is not a valid example number',ex_num)])
end
            


% Get the coefficients of f(x,y) in Bernstein form
fxy = GetCoefficients(sym_factor_arr_f);

% Get the symbolic polynomial in power basis
symbolic_fxy = GetSymbolicPoly(sym_factor_arr_f);

% Get the total degree of f(x,y)
m = double(feval(symengine, 'degree', symbolic_fxy));

fprintf('%s : Polynomial f(x,y). \n',mfilename);
display(symbolic_fxy);
end