function [root_mult_arr_f] = GetCoefficientsFromSymbolicRoots(root_mult_arr_f)

    sym_factor_arr_f = GetFactors(root_mult_arr_f);
    
    % Get the coefficients of f(x,y) in Bernstein form
    root_mult_arr_f = GetCoefficients(sym_factor_arr_f);

end