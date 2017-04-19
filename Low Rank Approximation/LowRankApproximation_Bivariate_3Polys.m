function [fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, alpha_lr, beta_lr, th1_lr, th2_lr] = LowRankApproximation_Bivariate_3Polys...
    (fxy, gxy, hxy, alpha, beta, th1, th2, k1, k2, idx_col)
% Compute low rank approximation of the Sylvester matrix S(f,g) either
% SNTLN or STLN.
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y) 
%
% hxy : (Matrix) Coefficients of polynomial h(x,y)
%
% alpha : (Float) \alpha
%
% beta : (Float) \beta
% 
% th1 : (Float) \theta_{1}
%
% th2 : (Float) \theta_{2}
%
% k1 : (Int) Degree of d(x,y) with respect to x
% 
% k2 : (Int) Degree of d(x,y) with respect to y
%
% idx_col : (Int) Index of optimal column to be removed from S(f,g)
%
% % Outputs
%
% fxy_lr : (Matrix) Coefficients of f(x,y) with perturbations
%
% gxy_lr : (Matrix) Coefficients of g(x,y) with perturbations
%
% alpha_lr : (Float) Refined \alpha
%
% beta_lr : (Float) Refined \beta
%
% th1_lr : (Float) Refined \theta_{1}
%
% th2_lr : (Float) Refined \theta_{2}


% Initialise global settings
global SETTINGS

switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    case 'Standard SNTLN'
        
        error([mfilename ' : Code Not completed'])
        
        % Apply SNTLN improvements
        %[ fxy_lr,gxy_lr,uxy_lr,vxy_lr,alpha_lr,th1_lr,th2_lr] = ...
         %   SNTLN(fxy, gxy, alpha, th1, th2, k1, k2, idx_col);
        
       
        
    case 'Standard STLN'

        error([mfilename ' : Code Not completed'])
        
        
        % Obtain polynomials in Modified Bernstein Basis, using initial values of
        % alpha and theta.

        % Multiply the rows of fxy_matrix by theta1, and multiply the cols of
        % fxy_matrix by theta2.
        fww = GetWithThetas(fxy, th1, th2);

        % Multiply the rows of gxy_matrix by theta1, and multiply the cols of
        % gxy_matrix by theta2.
        a_gww = alpha .* GetWithThetas(gxy, th1, th2);
        
        % Multiply 
        b_hww = beta .* GetWithTheats(hxy, th1, th2);
        
        % Perform STLN Computation.
        [fww_lr, a_gww_lr, uww_lr, vww_lr] = STLN_3Polys(fww, a_gww, b_hww, k1, k2,idx_col);
        
        % Scale outputs to obtain f(x,y) and g(x,y).
        fxy_lr = GetWithoutThetas(fww_lr,th1,th2);
        gxy_lr = GetWithoutThetas(a_gww_lr,th1,th2) ./ alpha;
        uxy_lr = GetWithoutThetas(uww_lr,th1,th2);
        vxy_lr = GetWithoutThetas(vww_lr,th1,th2);
        
        alpha_lr = alpha;
        beta_lr = beta;
        
        th1_lr = th1;
        th2_lr = th2;
        
        S1 = BuildDTQ_Bivariate_2Polys(fxy, gxy, k1, k2);
        S2 = BuildDTQ_bivar_2Polys(fxy_lr, gxy_lr, k1, k2);
        S3 = BuildDTQ_bivar_2Polys(fww, a_gww, k1, k2);
        S4 = BuildDTQ_bivar_2Polys(fww_lr, a_gww_lr, k1, k2);
        
        vSingularValues1 = svd(S1);
        vSingularValues2 = svd(S2);
        vSingularValues3 = svd(S3);
        vSingularValues4 = svd(S4);
        
       
        
    case 'None'
        % Dont Apply SNTLN improvements

        % Get polynomials in preprocessed form
        fww = GetWithThetas(fxy, th1, th2);
        a_gww = alpha .* GetWithThetas(gxy, th1, th2);
        b_hww = beta .* GetWithThetas(hxy, th1, th2);
        
        
        [uww, vww, www] = GetCofactors_Bivariate_3Polys(fww, a_gww, b_hww, k1, k2);
        
        % Get polynomials u(x,y), v(x,y) and w(x,y)
        uxy = GetWithoutThetas(uww, th1, th2);
        vxy = GetWithoutThetas(vww, th1, th2);
        wxy = GetWithoutThetas(www, th1, th2);
        
        
        uxy_lr = uxy;
        vxy_lr = vxy;
        wxy_lr = wxy;
        
        fxy_lr = fxy;
        gxy_lr = gxy;
        hxy_lr = hxy;
        
        alpha_lr = alpha;
        beta_lr = beta;
        th1_lr = th1;
        th2_lr = th2;
        
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
        
    otherwise
        error('bool_SNTLN is either (Standard SNTLN) or (None)')
end