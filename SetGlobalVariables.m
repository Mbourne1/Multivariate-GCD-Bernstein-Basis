function [] = SetGlobalVariables(mean_method,bool_alpha_theta,low_rank_approx_method)

% Set Variables

% BOOL_Q ('y'/'n')
%       'y' - Include the matrix Q in the Sylvester matrix such that it is
%       given by D^{-1}T(f,g)Q
%       'n' - Exclude the matrix Q from the Sylvester matrix such that it is
%       given by D^{-1}T(f,g)
global BOOL_Q
BOOL_Q = 'y';

global MEAN_METHOD
MEAN_METHOD = mean_method;

% bool_alpha_theta ('y'/'n')
%       'y' - Include the three preprocessing operations.
%       'n' - Exclude the three preprocessing operations.
global BOOL_ALPHA_THETA
BOOL_ALPHA_THETA = bool_alpha_theta;

% BOOL_NOISE ('y'/'n')
%       'y' - Add noise to the coefficients of input polynomial f
%       'n' - Use exact form of input polynomial f
global BOOL_NOISE
BOOL_NOISE = 'y';

% bool_SNTLN ('y'/'n')
%       'Standard SNTLN' - Include SNTLN
%       'None' - Exclude SNTLN
global LOW_RANK_APPROXIMATION_METHOD
LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

% bool_plotgraphs ('y'/'n')
%       y - plot graphs
%       n - exclude plotting
global PLOT_GRAPHS
PLOT_GRAPHS = 'y';

% seed - SEED Number for noise generation
global SEED
SEED = 1024;

%% Variables related to get degree ()
global THRESHOLD
THRESHOLD = 1;

%% Variables related to Low rank approximation
global MAX_ERROR_SNTLN
MAX_ERROR_SNTLN = 1e-14;

global MAX_ITERATIONS_SNTLN
MAX_ITERATIONS_SNTLN = 50;