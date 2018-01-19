function [] = Experiment3LowRankApproximation_2Polys(ex_num, bool_preproc)
% Experiment considers the alternate low rank approximation methods while
% other variables are kept constant.
%
% % Inputs
%
% ex_num : (String) Example Number
%
% bool_preproc : (Boolean) Whether preprocessing is included


close all;
clc;

% Upper and lower level of noise added to coefficients of f(x,y) and g(x,y)
emin = 1e-10;
emax = 1e-10;


apf_method = 'None';

% Factorisation Build Method
factorisation_build_method = 'HCG';

% Rank Revealing Metric
% 'Minimum Singular Values'
rank_revealing_metric = 'Minimum Singular Values';

% Sylvester Build Method
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_matrix_variant = 'DTQ';

% Set preprocessing related variables
switch bool_preproc
    case true
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
        
end


low_rank_approx_method_arr = {...
    'None', ...
    'Standard STLN'};

for i = 1 : 1 : length(low_rank_approx_method_arr)
    
    low_rank_approx_method = low_rank_approx_method_arr{i};
    
    o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, ...
        bool_alpha_theta, low_rank_approx_method, apf_method, ...
        sylvester_matrix_variant, factorisation_build_method, ...
        rank_revealing_metric)
end


end