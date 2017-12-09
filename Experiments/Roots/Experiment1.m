
close all; clc;

% Constants
ex_num = '3';
emin = 1e-12;
emax = 1e-10;
mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;
low_rank_approx_method = 'None';
apf_method = 'None';
sylvester_build_method = 'DTQ';
factorisation_build_method = 'HCG';
rank_revealing_metric = 'Minimum Singular Values';
nEquations = '2';

o_roots_Bivariate(ex_num, emin, emax, mean_method,...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_build_method, factorisation_build_method, rank_revealing_metric, nEquations)