function [] = Experiment2SylvesterVariants_3Polys(ex_num, bool_preproc)


close all; clc;

% Examples to not use
% 3 - This is univariate
% 4 - This is univariate
% 5 - Degree too low
% 6 - Degree of GCD too close to min(m1,n1,o1) min(m2,n2,o2)
% 8 - Degree t2 = min(m2,n2,k2)



emin = 1e-12;
emax = 1e-10;
low_rank_approx_method = 'None';
apf_method = 'None';

%ex_num_arr = {strcat(ex_num,'a'), strcat(ex_num,'b'), strcat(ex_num,'c')};
ex_num_arr = {strcat(ex_num,'a')};


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

%arrSylvesterFormats = {'T','DT','TQ', 'DTQ'};
arrSylvesterFormats = {'DTQ'};
nEquations = '2';

for i2 = 1 : 1 : length(ex_num_arr)
    
    ex_num = ex_num_arr{i2};
    
    for i = 1 : 1 : length(arrSylvesterFormats)
        
        
        sylvester_build_method = arrSylvesterFormats{i};
        
        o_gcd_Bivariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
            low_rank_approx_method, apf_method, sylvester_build_method, factorisation_build_method, ...
            rank_revealing_metric, nEquations)
        
    end
end


% nEquations = '3';
% o_gcd_Bivariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
%     low_rank_approx_method, apf_method, sylvester_build_method, factorisation_build_method, ...
%     rank_revealing_metric, nEquations)



end