function [] = o_roots_Bivariate(ex_num, emin, emax, mean_method,...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, factorisation_build_method, rank_revealing_metric, nEquations)
% o_roots_Bivariate(ex_num, el, mean_method, bool_alpha_theta, ...
%   low_rank_approx_method, apf_method, sylvester_matrix_variant, factorisation_build_method)
%
% Given an example number and set of parameters, obtain the roots of the
% example polynomial, where the polynomial is in the Bernstein form.
%
% Inputs
%
%
% ex_num : (String) Example Number
%
% el : (Sci) Lower noise level e.g 1e-10
%
% em : (Sci) Upper Noise level e.g 1e-12
%
% mean_method : (String)
%       'None'
%       'Geometric Mean Matlab Method'
%
% bool_alpha_theta (Boolean)
%       true - Include Preprocessing
%       false - Exclude Preprocessing
%
% low_rank_approx_method (String)
%       'Standard STLN' : Include STLN
%       'Standard SNTLN' : Include SNTLN
%       'None' - Exclude SNTLN
%
% apf_method : (String)
%       'None'
%       'Standard APF'
%
% sylvester_matrix_variant : (String)
%       'DTQ'
%       'DT'
%       'TQ'
%
% factorisation_build_method (String)
%       'HCG'
%       'HC'
%       'CG'
%
% rank_revealing_metric (String)
%
%
% % Examples
%
% >> o_roots_Bivariate('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard APF', 'DTQ', 'HCG')
% >> o_roots_Bivariate('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ', 'HCG')
% >> o_roots_Bivariate('1', 1e-10, 1e-12, 'None', true, 'None', 'None', 'DTQ', 'HCG')

% Clear path, and add all subfolders
addpath(genpath(pwd));

% Set the global variables
SetGlobalVariables_Roots(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
    factorisation_build_method, rank_revealing_metric, nEquations)

% Given the example number, return the coefficients of the bivariate
% polynomial f(x,y)
[fxy, arr_fxy_exact, arr_hxy_exact, arr_wxy_exact, ~] = ...
    Examples_Roots_Bivariate(ex_num);

% Add noise to the coefficients of polynomial f(x,y)
[fxy, ~] = AddVariableNoiseToPoly(fxy, emin, emax);

% Compute the roots of the polynomial
%
% 2 Poly GCD : Performs two sequences of polynomial GCD computations. The
% first sequence takes the GCD of a polynomial and its derivative with
% respect to x, then the second sequence takes the GCD of the polynomial 
% and its derivative with respect to y.
%
% 3 Poly GCD : Performs one sequence of polynomial GCD Computations. The
% set of GCD problems involve three polynomials; the polynomial f; its
% partial derivative with respect to x and its partial derivative with 
% respect to y
        
o_roots_mymethod_newmethod(fxy);
        



end