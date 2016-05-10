function [] = o_roots(ex_num,el,mean_method,bool_alpha_theta, low_rank_approx_method)
% o_roots(ex_num,el,mean_method,bool_alpha_theta, low_rank_approx_method)
%
% Given an example number and set of parameters, obtain the roots of the
% example polynomial, where the polynomial is in the Bernstein form.
%
% Inputs 
%
%
% ex_num : Example Number
%
% el     : Lower noise level
%
% mean_method : 
%       'None'
%       'Geometric Mean Matlab Method'
%
% bool_alpha_theta ('y'/'n')
%       y - Include Preprocessing
%       n - Exclude Preprocessing
%
% low_rank_approx_method ('y'/'n')
%       'Standard SNTLN' : Include SNTLN
%       'None' - Exclude SNTLN
%
% >> o_roots('1',1e-10,'Geometric Mean Matlab Method','y', 'None')

% Set the global variables
global SETTINGS

SetGlobalVariables(mean_method,bool_alpha_theta,low_rank_approx_method)

% %
%                   Get Example

% Given the example number, return the coefficients of the bivariate
% polynomial f(x,y)
[fxy_matrix] = Examples_Roots(ex_num);

% Plot the surface of the bivariate polynomial f(x,y)
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        PlotImplicitBezierSurface(fxy_matrix)
    case 'n'
        
    otherwise
        error('bool_plotgraphs must be either y or n')
end     

% Add noise to the coefficients
switch SETTINGS.BOOL_NOISE
    case 'y'
        % Add noise to the coefficients of polynomial f(x,y)
        [fxy_matrix,~] = Noise2(fxy_matrix,el);
    case 'n'
        % Dont add noise to coefficients of polynomial f(x,y)
    otherwise
        error('bool_noise must be either y or n')
end

o_roots_mymethod(fxy_matrix)


end