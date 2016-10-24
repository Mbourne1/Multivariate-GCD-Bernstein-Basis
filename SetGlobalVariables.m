function [] = SetGlobalVariables(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method)


global SETTINGS

%-------------------------------------------------------------------------
%
%   Other Settings
%

% Example Number
SETTINGS.EX_NUM = ex_num;

% bool_plotgraphs ('y'/'n')
%       y - plot graphs
%       n - exclude plotting

SETTINGS.PLOT_GRAPHS = 'y';


% ------------------------------------------------------------------------
%
%       Sylvester Matrix Settings
%
%

% BOOL_Q ('y'/'n')
%       'y' - Include the matrix Q in the Sylvester matrix such that it is
%       given by D^{-1}T(f,g)Q
%       'n' - Exclude the matrix Q from the Sylvester matrix such that it is
%       given by D^{-1}T(f,g)

SETTINGS.BOOL_Q = 'y';

% -------------------------------------------------------------------------
%
%       GCD Degree Computation Settings
%
%


%
% Total
% Relative
% Both
%
SETTINGS.CALC_METHOD = 'Relative';

% Variables related to get degree ()

SETTINGS.THRESHOLD = 1;

SETTINGS.THRESHOLD_RANK = 1;


%--------------------------------------------------------------------------
%
%       Preprocessing Settings
%

SETTINGS.MEAN_METHOD = mean_method;

% bool_alpha_theta ('y'/'n')
%       'y' - Include the three preprocessing operations.
%       'n' - Exclude the three preprocessing operations.

SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;


%-------------------------------------------------------------------------
%
%       Noise Settings
%
%


% BOOL_NOISE ('y'/'n')
%       'y' - Add noise to the coefficients of input polynomial f
%       'n' - Use exact form of input polynomial f
SETTINGS.BOOL_NOISE = 'y';


% seed - SEED Number for noise generation
SETTINGS.SEED = 1024;

% Set lower Noise level
SETTINGS.EMIN = emin;

%--------------------------------------------------------------------------
%
%       Low Rank Approx Settings
%
%
%

% LOW_RANK_APPROXIMATION_METHOD 
%       'Standard SNTLN' - Include SNTLN
%       'None' - Exclude SNTLN

SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;


% Variables related to Low rank approximation
SETTINGS.MAX_ERROR_SNTLN = 1e-14;

% Set maximum number of iterations for STLN
SETTINGS.MAX_ITERATIONS_SNTLN = 50;


end