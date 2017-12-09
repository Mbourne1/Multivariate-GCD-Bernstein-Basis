function [] = Experiment3LowRankApproximation_2Polys(ex_num, bool_preproc)
close all; clc;


emin = 1e-10;
emax = 1e-10;

apf_method = 'None';
factorisation_build_method = 'HCG';
rank_revealing_metric = 'Minimum Singular Values';

switch bool_preproc
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
    case false
        mean_method = 'None';
        bool_alpha_theta = false;
end

sylvester_build_method = 'DTQ';



low_rank_approx_method = 'None';
o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, factorisation_build_method, rank_revealing_metric)

low_rank_approx_method = 'Standard STLN';
o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, factorisation_build_method, rank_revealing_metric)

%low_rank_approx_method = 'Standard SNTLN';
%o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, factorisation_build_method, rank_revealing_metric)
end