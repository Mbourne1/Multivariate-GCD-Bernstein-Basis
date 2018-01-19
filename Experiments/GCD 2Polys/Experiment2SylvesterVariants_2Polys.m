function [] = Experiment2SylvesterVariants_2Polys(ex_num, bool_preproc)
% Compute the GCD using different variants of the Sylvester subresultant
% matrices.
%
% % Inputs
%
% ex_num : (String) Example number
%
% bool_preproc : (Boolean) 
%
%
% % Examples
%
% Experiment2SylvesterVariants_2Polys('1');


close all; clc;

% Set noise levels
emin = 1e-9;
emax = 1e-9;

% Low Rank Approximation Method
low_rank_approx_method = 'None';

% Approximate Polynomial Factorisation Method
apf_method = 'None';

% Factorisation build method
factorisation_build_method = 'HCG';

% Rank Revealing Metric
% 'Minimum Singular values'
rank_revealing_metric = 'Minimum Singular Values';


switch bool_preproc
    
    case true
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
        
end

% Get array of subresutlant matrix variants
arrSylvesterSubresultantVariants = {'T', ...
    'DT', ...
    'TQ', ...
    'DTQ'};

degree_method = 'All Subresultants';

% For each subresultant matrix variant
for i = 1 : 1 : length(arrSylvesterSubresultantVariants)
    
    % Set variant
    sylvester_matrix_variant = arrSylvesterSubresultantVariants{i};
    
    % Get GCD
    o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, ...
        bool_alpha_theta, low_rank_approx_method, apf_method, ...
        sylvester_matrix_variant, factorisation_build_method, ...
        rank_revealing_metric, degree_method)
    
    
end



end