function [] = Experiment2SylvesterVariants_2Polys(ex_num)
% Compute the GCD using different variants of the Sylvester subresultant
% matrices.
%
% % Inputs
%
% ex_num : (String) Example number
%
%
% % Examples
%
% Experiment2SylvesterVariants_2Polys('1');


% Notes

%
%
% 1 : Too small%
% 2 : All methods pass, too easy
% 3 : too small
% 4 : too small
% 5 : too small
% 6 :
%


close all;
clc;

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


bool_preproc = true;
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

% For each subresultant matrix variant
for i = 1 : 1 : length(arrSylvesterSubresultantVariants)
    
    % Set variant
    sylvester_build_method = arrSylvesterSubresultantVariants{i};
    
    % Get GCD
    o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, ...
        bool_alpha_theta, low_rank_approx_method, apf_method, ...
        sylvester_build_method, factorisation_build_method, ...
        rank_revealing_metric)
    
    
end



end