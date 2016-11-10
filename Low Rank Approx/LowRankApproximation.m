function [fxy_matrix,gxy_matrix,alpha,th1,th2] = LowRankApproximation...
    (fxy_matrix_n,gxy_matrix_n,alpha,th1,th2,k1,k2,idx_col)
% Compute low rank approximation of the Sylvester matrix S(f,g) either
% SNTLN or STLN.
%
% % Inputs.
%
% fxy_matrix_n : Coefficients of polynomial f(x,y)
%
% gxy_matrix_n : Coefficients of polynomial g(x,y)
%
% alpha : \alpha
%
% th1 : \theta_{1}
%
% th2 : \theta_{2}
%
% k1 : Degree of d(x,y) with respect to x
% 
% k2 : Degree of d(x,y) with respect to y
%
% idx_col : Index of optimal column to be removed from S(f,g)
%
% % Outputs
%
% fxy_matrix : Coefficients of f(x,y) with perturbations
%
% gxy_matrix : Coefficients of g(x,y) with perturbations
%
% alpha : Refined \alpha
%
% th1 : Refined \theta_{1}
%
% th2 : Refined \theta_{2}


% Initialise global settings
global SETTINGS

switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    case 'Standard SNTLN'
        
        
        % Apply SNTLN improvements
        [ fxy_matrix,gxy_matrix,alpha,th1,th2,x] = ...
            SNTLN( fxy_matrix_n,gxy_matrix_n, alpha, th1, th2,k1,k2, idx_col);
        
        fprintf('Input Polynomial f(x,y)')
        
    case 'Standard STLN'

        % Obtain polynomials in Modified Bernstein Basis, using initial values of
        % alpha and theta.

        % Multiply the rows of fxy_matrix by theta1, and multiply the cols of
        % fxy_matrix by theta2.
        fww_matrix = GetWithThetas(fxy_matrix_n,th1,th2);

        % Multiply the rows of gxy_matrix by theta1, and multiply the cols of
        % gxy_matrix by theta2.
        gww_matrix = GetWithThetas(gxy_matrix_n,th1,th2);
        
        % Perform STLN Computation.
        [fww_matrix, gww_matrix, ~] = STLN(fww_matrix, alpha.*gww_matrix,k1,k2,idx_col);
        
        % Scale outputs to obtain f(x,y) and g(x,y).
        fxy_matrix = GetWithoutThetas(fww_matrix,th1,th2);
        gxy_matrix = GetWithoutThetas(gww_matrix,th1,th2) ./ alpha;
        
        
        
    case 'None'
        % Dont Apply SNTLN improvements
        fxy_matrix = fxy_matrix_n;
        gxy_matrix = gxy_matrix_n;

        
    otherwise
        error('bool_SNTLN is either (Standard SNTLN) or (None)')
end