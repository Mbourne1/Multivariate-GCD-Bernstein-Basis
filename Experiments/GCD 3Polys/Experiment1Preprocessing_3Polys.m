function [] = Experiment1Preprocessing_3Polys(ex_num, bool_preproc)
%
% % Inputs
%
% ex_num : (String) Example number
%
% bool_preproc : (Boolean) Determine whether to preprocess the polynomials
% f(x,y), g(x,y) and h(x,y) in the computation of the GCD.
%
% % Examples
%
% >> Experiment1Preprocessing_3Polys('1', true);
%
% >> Experiment1Preprocessing_3Polys('1', false);


% 3 - Bad Example - Univariate
% 4 - Bad Example - Univariate
% 5 - Small Example

% Set upper and lower bound for noise added to coefficients of f(x,y) and
% g(x,y)
emin = 1e-10;
emax = 1e-8;

% Low Rank Approximation Method
low_rank_approx_method = 'None';


apf_method = 'None';

% Sylvester Subresultant Matrix Variants
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_matrix_variant = 'DTQ';


factorisation_build_method = 'HCG';

% Rank Revealing Metric
% 'Minimum Singular Values'
rank_revealing_metric = 'Minimum Singular Values';

% nEquations determines the structure of the subresultant matrices. 
% 2 : Gives a (2x3) partitioned subresultant matrix
% 3 : Gives a (3x3) partitioned subresultant matrix
nEquations = '3';


% % Degree Method
% 'Total'
% 'All Subresultants'
% 'Linear Method'

degree_method = 'All Subresultants';

% Set preprocessing related variables
switch bool_preproc
    case true
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
        
end

o_gcd_Bivariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_matrix_variant, ...
    factorisation_build_method, rank_revealing_metric, nEquations, ...
    degree_method)

end
