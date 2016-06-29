function [fxy_matrix,gxy_matrix,alpha,th1,th2] = LowRankApproximation...
    (fxy_matrix_n,gxy_matrix_n,alpha,th1,th2,t1,t2,opt_col)


global SETTINGS

switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    case 'Standard SNTLN'
        
        
        % Apply SNTLN improvements
        [ fxy_matrix,gxy_matrix,alpha,th1,th2,x] = ...
            SNTLN( fxy_matrix_n,gxy_matrix_n, alpha, th1, th2,t1,t2, opt_col);
        
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
        [fww_matrix, gww_matrix, ~] = STLN(fww_matrix, alpha.*gww_matrix,t1,t2,opt_col);
        
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