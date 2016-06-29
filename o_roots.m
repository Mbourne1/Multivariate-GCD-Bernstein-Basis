function [] = o_roots(ex_num,emin,emax,mean_method,bool_alpha_theta, low_rank_approx_method)
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
% em : Upper Noise level 
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
% >> o_roots('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', 'y', 'None')

% Set the global variables
global SETTINGS

SetGlobalVariables(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method)

% %
%                   Get Example

% Given the example number, return the coefficients of the bivariate
% polynomial f(x,y)
[fxy_matrix,M] = Examples_Roots(ex_num);

% Plot the surface of the bivariate polynomial f(x,y)
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        try
        PlotImplicitBezierSurface(fxy_matrix)
        catch
        end
    case 'n'
        
    otherwise
        error('bool_plotgraphs must be either y or n')
end     


% Add noise to the coefficients of polynomial f(x,y)
[fxy_matrix,~] = Noise2(fxy_matrix,emin,emax);


o_roots_mymethod(fxy_matrix,M)


end