function [] = SetGlobalVariables_GCD_3Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, factorisation_build_method, rank_revealing_metric, nEquations )
%
% % Inputs
%
% ex_num : (String)
%
% emin : (Float)
%
% emax : (Float)
%
% mean_method : (String)
%
% bool_alpha_theta : (Boolean)
%
% low_rank_approx_method : (String) 
%
% apf_method : (String)
%
% sylvester_matrix_variant : (String) 
%
% factorisation_build_method : (String)
%
% rank_revealing_metric : (String)
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Minimum Singular Values
%   * Residuals

global SETTINGS

%-------------------------------------------------------------------------
%
%   Other Settings
%

% Example Number
SETTINGS.EX_NUM = ex_num;

% bool_plotgraphs (Boolean)
%   * true : Plot graphs
%   * false : Exclude plotting
SETTINGS.PLOT_GRAPHS = true;


% ------------------------------------------------------------------------
%
% Sylvester Matrix Settings
%
%

SETTINGS.SYLVESTER_MATRIX_VARIANT = sylvester_matrix_variant;
SETTINGS.SYLVESTER_MATRIX_3POLY_N_EQUATIONS = nEquations;

% Version 1 - Read antidiagonals
% Version 2 - Read down columns from left to right
SETTINGS.VECTORISATION_METHOD = 'Version 1';
%-------------------------------------------------------------------------
%
% Factorisation Matrix Settings
%

SETTINGS.FACTORISATION_BUILD_METHOD = factorisation_build_method;

% -------------------------------------------------------------------------
%
% GCD Degree Computation Settings
%
%

% Variables related to get degree ()

SETTINGS.THRESHOLD = 1;

SETTINGS.THRESHOLD_RANK = 1;

% Metric used to compute the degree of the GCD
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Minimum Singular Values
%   * Residuals
SETTINGS.RANK_REVEALING_METRIC = rank_revealing_metric;

%--------------------------------------------------------------------------
%
%       Preprocessing Settings
%

SETTINGS.MEAN_METHOD = mean_method;

% bool_alpha_theta true/false
%   * true :  Include the three preprocessing operations.
%   * false :  Exclude the three preprocessing operations.

SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;


%-------------------------------------------------------------------------
%
% Noise Settings
%
%
% seed - SEED Number for noise generation
SETTINGS.SEED = 1024;

% Set lower Noise level
SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;
%--------------------------------------------------------------------------
%
% Low Rank Approx Settings
%

% LOW_RANK_APPROXIMATION_METHOD 
%       'Standard SNTLN' - Include SNTLN
%       'None' - Exclude SNTLN

SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;


% Variables related to Low rank approximation
SETTINGS.MAX_ERROR_SNTLN = 1e-11;

% Set maximum number of iterations for STLN
SETTINGS.MAX_ITERATIONS_SNTLN = 50;

%-------------------------------------------------------------------------
% 
% APF METHOD

SETTINGS.APF_METHOD = apf_method;

%-------------------------------------------------------------------------

end