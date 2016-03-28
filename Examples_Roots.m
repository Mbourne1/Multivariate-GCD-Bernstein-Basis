function [fxy_matrix_exact,m] = Examples_Roots(ex_num)
% Examples_Roots(ex_num)
%
% Given an example number, return the matrix of coefficients of a
% polynomial f(x,y)
%
% Inputs.
%
% ex_num : Example Number
%
%



root_mult_array_f_x = [];
root_mult_array_f_y = [];
factor_f_xy = [];

switch ex_num
    case '1'
        % Relates to roots_examples tex file example 1
        
        % Get roots of f with respect to x
        root_mult_array_f_x =...
            [
            0.9   2;
            ];
        
        % Get roots of f with respect to y
        root_mult_array_f_y = ...
            [
            0.3   1
            ];
        
        % Get roots of f with respect to x and y
        % (x+y-0.1)
        factor_f_xy{1,1} = ...
            [
            0.1  1.1;
            1.1  2.1
            ];
        
        m = 5;
        
        
        
    case '2'
        
        % Get roots and multiplicities wrt x
        root_mult_array_f_x = [...
            0.8     1;
            0.1     2;
            0.5     2;
            0.9     3;
            1.1234  4;
            ];
        
        m = 12;
        
    case '3'
        % Roots with respect to x
        root_mult_array_f_x = ...
            [
            0.5 1
            0.1 2
            ];
        % Roots with respect to y
        root_mult_array_f_y = ...
            [
            0.4 2
            ];
        
        m = 5;
        
        
    case '4'
        % (x-0.5)(x-0.1)^2 (x-0.7972)^3 (y-0.4)^2 (y-0.687) (x+y-0.1)
        
        % Roots with respect to x
        root_mult_array_f_x = ...
            [
            0.5     1
            0.1     2
            0.7972  3
            ];
        % Roots with respect to x
        root_mult_array_f_y = ...
            [
            0.4     2
            0.687   1
            ];
        % Roots with respect to x and y
        factor_f_xy{1,1} = ...
            [
            0.1  1.1;
            1.1  2.1
            ];
        
        factor_f_xy{2,1} = ...
            [
            0.1  1.1;
            1.1  2.1
            ];
        
        
        m = 13;
    case '5'
        % (x-0.5)(x-0.1)^2 (x-0.7972)^3 (y-0.4)^2 (y-0.687) (x+y-0.1)
        
        % Roots with respect to x
        root_mult_array_f_x = ...
            [
            0.5     1
            0.1     2
            0.7972  3
            ];
        % Roots with respect to y
        root_mult_array_f_y = ...
            [
            0.4     2
            0.687   1
            ];
        % Roots with respect to x and y
        factor_f_xy{1,1} = ...
            [
            0.1  1.1;
            1.1  2.1
            ];
        
        factor_f_xy{2,1} = ...
            [
            0.1  1.1;
            1.1  2.1
            ];
        factor_f_xy{3,1} = ...
            [
            0.1  1.1;
            1.1  2.1
            ];
        
        m = 15;
        
    case '6'
        
        % Roots with respect to x
        root_mult_array_f_x = ...
            [
            0.5     1
            0.1     2
            0.7972  3
            0.2475  4
            ];
        % Roots with respect to y
        root_mult_array_f_y = ...
            [
            0.4     2
            0.687   1
            ];
        
        % Roots with respect to x and y
        factor_f_xy{1,1} = ...
            [
            0.1  1.1;
            1.1  2.1
            ];
        
        factor_f_xy{2,1} = ...
            [
            0.1  1.1;
            1.1  2.1
            ];
        factor_f_xy{3,1} = ...
            [
            0.1  1.1;
            1.1  2.1
            ];
        
        m = 16;
        
    otherwise
        error('err')
end

% Print out the roots
PrintRoots_Bivariate(root_mult_array_f_x,root_mult_array_f_y,factor_f_xy);

% Get the array of factors of f (in x)
factors_x = mult_roots_x(root_mult_array_f_x);

% Get the array of factors of f (in y)
factors_y = mult_roots_y(root_mult_array_f_y);

% Collect all factors
factors_f = [factors_x; factors_y; factor_f_xy];

% Get the coefficients of the polynomail f(x,y) from the set of roots
fxy_matrix_exact = BuildPoly_NonSeparable(factors_f);

end





