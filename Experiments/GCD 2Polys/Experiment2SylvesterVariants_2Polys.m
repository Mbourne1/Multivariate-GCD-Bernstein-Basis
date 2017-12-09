function [] = Experiment2SylvesterVariants_2Polys(ex_num)
% Compute the 
%
%
% % Inputs
%
% ex_num : (String) 
%

%
%
% 1 : Too small% 
% 2 : All methods pass, too easy
% 3 : too small
% 4 : too small
% 5 : too small
% 6 : 
%
%
%
%
%
%




close all; clc;

emin = 1e-9;
emax = 1e-9;
low_rank_approx_method = 'None';
apf_method = 'None';

factorisation_build_method = 'HCG';
rank_revealing_metric = 'Minimum Singular Values';
%mean_method = 'Geometric Mean Matlab Method';
%bool_alpha_theta = true;

mean_method = 'None';
bool_alpha_theta = false;


arrSylvesterFormats = {'T', 'DT', 'TQ', 'DTQ'};


for i = 1 : 1 : length(arrSylvesterFormats)

    
    sylvester_build_method = arrSylvesterFormats{i};
    
    o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, factorisation_build_method, rank_revealing_metric)
    
    
end



end