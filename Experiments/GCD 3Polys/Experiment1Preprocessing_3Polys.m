function [] = Experiment1Preprocessing_3Polys(ex_num, bool_preproc)

%
% Examples
%
% 3 - Bad Example - Univariate
% 4 - Bad Example - Univariate
% 5 - Small Example
% 6 - 
%

%ex_num = '13a';
emin = 1e-10;
emax = 1e-8;
low_rank_approx_method = 'None';
apf_method = 'None';
sylvester_build_method = 'DTQ';
factorisation_build_method = 'HCG';
rank_revealing_metric = 'Minimum Singular Values';

nEquations = '3';

%degree_method = 'Total';
degree_method = 'All Subresultants';
%degree_method = 'Linear Method';

switch bool_preproc
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
    case false
        mean_method = 'None';
        bool_alpha_theta = false;
end

o_gcd_Bivariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, factorisation_build_method, rank_revealing_metric, nEquations, degree_method)

end
