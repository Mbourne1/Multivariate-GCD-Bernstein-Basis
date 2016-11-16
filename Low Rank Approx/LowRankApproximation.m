function [fxy_lr, gxy_lr, uxy_lr,vxy_lr, alpha_lr, th1_lr, th2_lr] = LowRankApproximation...
    (fxy, gxy, alpha, th1, th2, m, n, k, k1, k2, idx_col)
% Compute low rank approximation of the Sylvester matrix S(f,g) either
% SNTLN or STLN.
%
% % Inputs.
%
% fxy : Coefficients of polynomial f(x,y) in the Bernstein basis
%
% gxy : Coefficients of polynomial g(x,y) in the Bernstein basis
%
% alpha : \alpha
%
% th1 : \theta_{1}
%
% th2 : \theta_{2}
%
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomial g(x,y)
%
% k : Total degree of polynomial d(x,y)
%
% k1 : Degree of d(x,y) with respect to x
% 
% k2 : Degree of d(x,y) with respect to y
%
% idx_col : Index of optimal column to be removed from S(f,g)
%
% % Outputs
%
% fxy_lr : Coefficients of f(x,y) with perturbations
%
% gxy_lr : Coefficients of g(x,y) with perturbations
%
% alpha_lr : Refined \alpha
%
% th1_lr : Refined \theta_{1}
%
% th2_lr : Refined \theta_{2}


% Initialise global settings
global SETTINGS

switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    case 'Standard SNTLN'
        
        
        % Apply SNTLN improvements
        [ fxy_lr,gxy_lr,uxy_lr,vxy_lr,alpha_lr,th1_lr,th2_lr] = ...
            SNTLN(fxy, gxy, alpha, th1, th2, k1, k2, idx_col);
        
       
        
    case 'Standard STLN'

        % Obtain polynomials in Modified Bernstein Basis, using initial values of
        % alpha and theta.

        % Multiply the rows of fxy_matrix by theta1, and multiply the cols of
        % fxy_matrix by theta2.
        fww = GetWithThetas(fxy,th1,th2);

        % Multiply the rows of gxy_matrix by theta1, and multiply the cols of
        % gxy_matrix by theta2.
        a_gww = alpha .* GetWithThetas(gxy,th1,th2);
        
        % Perform STLN Computation.
        [fww_lr, a_gww_lr, uww_lr, vww_lr] = STLN(fww, a_gww,k1,k2,idx_col);
        
        % Scale outputs to obtain f(x,y) and g(x,y).
        fxy_lr = GetWithoutThetas(fww_lr,th1,th2);
        gxy_lr = GetWithoutThetas(a_gww_lr,th1,th2) ./ alpha;
        uxy_lr = GetWithoutThetas(uww_lr,th1,th2);
        vxy_lr = GetWithoutThetas(vww_lr,th1,th2);
        
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        S1 = BuildDTQ(fxy,gxy,k1,k2);
        S2 = BuildDTQ(fxy_lr,gxy_lr,k1,k2);
        S3 = BuildDTQ(fww,a_gww,k1,k2);
        S4 = BuildDTQ(fww_lr,a_gww_lr,k1,k2);
        
        vSingularValues1 = svd(S1);
        vSingularValues2 = svd(S2);
        vSingularValues3 = svd(S3);
        vSingularValues4 = svd(S4);
        
        figure()
        plot(log10(vSingularValues1),'-s','DisplayName','fxy,gxy');
        hold on
        plot(log10(vSingularValues2),'-s','DisplayName','fxy_lr,gxy_lr');
        plot(log10(vSingularValues3),'-s','DisplayName','fww gww');
        plot(log10(vSingularValues4),'-s','DisplayName','fww_lr gww_lr');
        hold off
        
    case 'None'
        % Dont Apply SNTLN improvements

        [uxy,vxy] = GetCofactors(fxy,gxy,m,n,k,k1,k2,alpha,th1,th2);
        
        uxy_lr = uxy;
        vxy_lr = vxy;
        fxy_lr = fxy;
        gxy_lr = gxy;
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        
    otherwise
        error('bool_SNTLN is either (Standard SNTLN) or (None)')
end