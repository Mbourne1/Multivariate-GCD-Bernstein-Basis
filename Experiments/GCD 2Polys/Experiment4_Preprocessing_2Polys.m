function [] = Experiment4_Preprocessing_2Polys(ex_num)
% This example considers the two different methods for the computation of
% the degree of the GCD of two bivariate polynomials in Bernstein form.
%
%
%
% % Inputs
%
% ex_num : (String) Example number
%
% bool_Preproc : (Boolean) Determine whether polynomials f(x,y) and g(x,y)
% are preprocessed in the computation of their GCD.
%
% % Examples
%
% Experiment1Preprocessing_2Polys('2', true)


close all;
clc;


% Examples

% 1
% 2 - Good Example - !!!(Used In Thesis)!!!
% 3 - Bad Example - univariate
% 4 - Bad Example - univariate
% 5 - Bad Example - too small
% 6 - Bad Example - t1 or t2 too close to min(m1, n1) min(m2, n2)
% 7 -
% 8 - Bad Example - t1 or t2 too close
% 9 - Good Example - 1e-8, 1e-6
% 10 - Bad Example -Univariate Example
% 11 - Good Example - 1e-8, 1e-6
% 12 - Good Example - 1e-8, 1e-6
% 13 - Good Example - 1e-9, 1e-9
% 14 - Great Example - 1e-10, 1e-8
% 15 - Great Example
% 16 - Good Example
% 17 -
% 18 -
% 19 -
% 20 -
% 21 - Good Example - 1e-10, 1e-10
% 22 -

% Good Examples

% ex_num = '2'; emin = 1e-8; emax = 1e-6;

% ex_num = '12'; emin = 1e-8; emax = 1e-10;
% ex_num = '21'; emin = 1e-10; emax = 1e-10;





% Set max and min noise levels
emin = 1e-6;
emax = 1e-6;

% Low Rank Approximation Method
% 'None'
low_rank_approx_method = 'None';

% Approximate Polynomial Factorisation Method
apf_method = 'None';

% Sylvester Build Method
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
sylvester_build_method = 'DTQ';

% Factorisation Build Method
% 'HCG'
factorisation_build_method = 'HCG';

% Rank Revealing Metric
% 'Minimum Singular Values'
rank_revealing_metric = 'Minimum Singular Values';





bool_preproc_arr = {true, false};


% Perform GCD computation with bool_preproc set to true then false
for i = 1 : 1 : length(bool_preproc_arr)
    
    % Set bool_preproc true or false
    bool_preproc = bool_preproc_arr{i};
    
    
    % Set preprocessing related variables
    switch bool_preproc
        case true
            
            mean_method = 'Geometric Mean Matlab Method';
            bool_alpha_theta = true;
            %bool_alpha_theta = false;
            
        case false
            
            mean_method = 'None';
            bool_alpha_theta = false;
            
        otherwise
            error("Error")
    end
    
    
    
    try
        tic;
        degree_method = 'Total';
        o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, ...
            bool_alpha_theta, low_rank_approx_method, apf_method, ...
            sylvester_build_method, factorisation_build_method, ...
            rank_revealing_metric, degree_method)
        toc;
    catch
        fprintf('Error')
    end
    
    
    
    
    tic;
    degree_method = 'All Subresultants';
    o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta,...
        low_rank_approx_method, apf_method, sylvester_build_method, ...
        factorisation_build_method, rank_revealing_metric, degree_method)
    toc
    
    
end






end