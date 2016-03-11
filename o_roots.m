function [] = o_roots(ex_num,el,bool_preproc, bool_sntln,seed)
%
% Given an example number and set of parameters, obtain the roots of the
% example polynomial, where the polynomial is in the Bernstein form.
% %                             Inputs 
% ex_num - Example Number
%
% el - Lower noise level
%
% BOOL_PREPROC ('y'/'n')
%       y - Include Preprocessing
%       n - Exclude Preprocessing
%
% BOOL_SNTLN ('y'/'n')
%       y - Include SNTLN
%       n - Exclude SNTLN
%
% SEED - SEED Number for noise generation
%
%
%%
%                           Set Variables

% bool_Q ('y'/'n')
%       y - Include the matrix Q in the Sylvester matrix such that it is
%       given by D^{-1}T(f,g)Q
%       n - Exclude the matrix Q from the Sylvester matrix such that it is
%       given by D^{-1}T(f,g)
global bool_Q
bool_Q = 'y';

% bool_preproc ('y'/'n')
%       y - Include the three preprocessing operations.
%       n - Exclude the three preprocessing operations.
global BOOL_PREPROC
BOOL_PREPROC = bool_preproc;

% bool_noise ('y'/'n')
%       y - Add noise to the coefficients of input polynomial f
%       n - Use exact form of input polynomial f
global BOOL_NOISE
BOOL_NOISE = 'y';

% bool_SNTLN ('y'/'n')
%       y - Include SNTLN
%       n - Exclude SNTLN
global BOOL_SNTLN
BOOL_SNTLN = bool_sntln;

% bool_plotgraphs ('y'/'n')
%       y - plot graphs
%       n - exclude plotting
global PLOT_GRAPHS
PLOT_GRAPHS = 'y';

% seed - SEED Number for noise generation
global SEED
SEED = seed;

%%

%                   Get Example

% Given the example number, return the coefficients of the bivariate
% polynomial f(x,y)
[fxy_matrix] = Examples_Roots(ex_num);

% Get dimensions of polynomial f(x,y)
[r,c] = size(fxy_matrix);
m1 = r-1;
m2 = c-1;

% Plot the surface of the bivariate polynomial f(x,y)
switch PLOT_GRAPHS
    case 'y'
        PlotImplicitBezierSurface(fxy_matrix)
    case 'n'
        
    otherwise
        error('bool_plotgraphs must be either y or n')
end     

% Add noise to the coefficients
switch BOOL_NOISE
    case 'y'
        % Add noise to the coefficients of polynomial f(x,y)
        [fxy_matrix,~] = Noise2(fxy_matrix,el);
    case 'n'
        % Dont add noise to coefficients of polynomial f(x,y)
    otherwise
        error('bool_noise must be either y or n')
end

GetRoots(fxy_matrix)


end