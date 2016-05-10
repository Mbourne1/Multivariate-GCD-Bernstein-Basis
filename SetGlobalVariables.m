function [] = SetGlobalVariables(mean_method,bool_alpha_theta,low_rank_approx_method)

% Set Variables

% BOOL_Q ('y'/'n')
%       'y' - Include the matrix Q in the Sylvester matrix such that it is
%       given by D^{-1}T(f,g)Q
%       'n' - Exclude the matrix Q from the Sylvester matrix such that it is
%       given by D^{-1}T(f,g)
global SETTINGS
SETTINGS.BOOL_Q = 'y';


SETTINGS.MEAN_METHOD = mean_method;

% bool_alpha_theta ('y'/'n')
%       'y' - Include the three preprocessing operations.
%       'n' - Exclude the three preprocessing operations.

SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;

% BOOL_NOISE ('y'/'n')
%       'y' - Add noise to the coefficients of input polynomial f
%       'n' - Use exact form of input polynomial f

SETTINGS.BOOL_NOISE = 'y';

% bool_SNTLN ('y'/'n')
%       'Standard SNTLN' - Include SNTLN
%       'None' - Exclude SNTLN

SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

% bool_plotgraphs ('y'/'n')
%       y - plot graphs
%       n - exclude plotting

SETTINGS.PLOT_GRAPHS = 'y';

% seed - SEED Number for noise generation
SETTINGS.SEED = 1024;

%% Variables related to get degree ()
SETTINGS.THRESHOLD = 1;

%% Variables related to Low rank approximation
SETTINGS.MAX_ERROR_SNTLN = 1e-14;
SETTINGS.MAX_ITERATIONS_SNTLN = 50;