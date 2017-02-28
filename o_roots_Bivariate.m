function [] = o_roots_Bivariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method)
% o_roots(ex_num,el,mean_method,bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method)
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
% bool_alpha_theta (true/false)
%       true - Include Preprocessing
%       false - Exclude Preprocessing
%
% low_rank_approx_method ('y'/'n')
%       'Standard STLN' : Include STLN
%       'Standard SNTLN' : Include SNTLN
%       'None' - Exclude SNTLN
%
% apf_method
%       'None'
%       'Standard APF'
%
% sylvester_build_method :
%       'DTQ'
%       'DT'
%       'TQ'
%
% >> o_roots('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', 'y', 'Standard STLN','Standard APF','DTQ')
% >> o_roots('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', 'y', 'None','None','DTQ')
% >> o_roots('1', 1e-10, 1e-12, 'None', 'n', 'None','None','DTQ')

% Add subfolders
restoredefaultpath

addpath(...
    'APF',...
    'Basis Conversion',...
    'Bernstein Methods',...
    'Build Matrices',...
    'Deconvolution',...
    'Examples',...
    'Examples/Examples Roots',...
    'Formatting',...
    'Get Cofactors',...
    'Get GCD Coefficients',...
    'Get GCD Degree',...
    'Low Rank Approx',...
    'Plotting',...
    'Preprocessing',...
    'Results',...
    'Root Finding Methods',...
    'Sylvester Matrix');


% Set the global variables
global SETTINGS
SetGlobalVariables(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method)


% Given the example number, return the coefficients of the bivariate
% polynomial f(x,y)
[fxy_matrix,M] = Examples_Roots(ex_num);

% Plot the surface of the bivariate polynomial f(x,y)
if( SETTINGS.PLOT_GRAPHS)
    
        try
            %PlotImplicitBezierSurface(fxy_matrix)
        catch
        end
    
end


% Add noise to the coefficients of polynomial f(x,y)
[fxy_matrix,~] = AddVariableNoiseToPoly(fxy_matrix,emin,emax);

switch root_finding_method
    case '2 Poly GCD'
        
        o_roots_mymethod(fxy_matrix,M);
        
    case '3 Poly GCD'
        
        o_roots_mymethod_newmethod(fxy_matrix,M);
        
    otherwise
        
        error([mfilename ': Error \n'])
end


end