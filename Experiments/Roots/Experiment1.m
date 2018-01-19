function [] = Experiment1(ex_num, bool_preproc)
% 
% % Inputs
% 
% ex_num : (String) Example number
%
% bool_preproc : (Boolean) Determine whether subresultant matrices are
% preprocessed in each GCD computation.

close all;
clc;

% Constants

% Set upper and lower limit of noise added to the coefficients of f(x,y)
emin = 1e-12;
emax = 1e-10;

% Set preprocessing related variables
switch bool_preproc
    case true
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
end



low_rank_approx_method = 'None';

apf_method = 'None';

% Sylvester Subresultant Matrix Variants
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_matrix_variant = 'DTQ';

factorisation_build_method = 'HCG';

rank_revealing_metric = 'Minimum Singular Values';

nEquations = '2';

o_roots_Bivariate(ex_num, emin, emax, mean_method,...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, factorisation_build_method, rank_revealing_metric, nEquations)



end