function [fxy_lr, gxy_lr, uxy_lr,vxy_lr, alpha_lr, th1_lr, th2_lr] = LowRankApproximation_Bivariate_2Polys...
    (fxy, gxy, alpha, th1, th2, k1, k2, idx_col)
% Compute low rank approximation of the Sylvester matrix S(f,g) either
% SNTLN or STLN.
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) in the Bernstein basis
%
% gxy : (Matrix) Coefficients of polynomial g(x,y) in the Bernstein basis
%
% alpha : (Float) \alpha
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
% th1_lr : (Float) Refined \theta_{1}
%
% th2_lr : (Float) Refined \theta_{2}


% Initialise global settings
global SETTINGS

switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    case 'Standard SNTLN'
        
        
        % Apply SNTLN improvements
        [ fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
            SNTLN(fxy, gxy, alpha, th1, th2, k1, k2, idx_col);
        
        
        fww = GetWithThetas(fxy, th1, th2);
        a_gww = alpha.*GetWithThetas(gxy, th1, th2);
        fww_lr = GetWithThetas(fxy_lr, th1_lr, th2_lr);
        a_gww_lr = alpha_lr .* GetWithThetas(gxy_lr, th1_lr, th2_lr);
        
        if SETTINGS.PLOT_GRAPHS_LRA == true
        PlotSingularValues(fxy, gxy, fxy_lr, gxy_lr, fww, a_gww, fww_lr, a_gww_lr)
        end
        
        
    case 'Standard STLN'
        
        % Obtain polynomials in Modified Bernstein Basis, using initial values of
        % alpha and theta.
        
        % Multiply the rows of fxy_matrix by theta1, and multiply the cols of
        % fxy_matrix by theta2.
        fww = GetWithThetas(fxy, th1, th2);
        a_gww = alpha .* GetWithThetas(gxy, th1, th2);
        
        % Perform STLN Computation.
        [fww_lr, a_gww_lr, uww_lr, vww_lr] = STLN(fww, a_gww, k1, k2, idx_col);
        
        % Scale outputs to obtain f(x,y) and g(x,y).
        fxy_lr = GetWithoutThetas(fww_lr, th1,th2);
        gxy_lr = GetWithoutThetas(a_gww_lr, th1, th2) ./ alpha;
        
        % Get u(x,y) and v(x,y)
        uxy_lr = GetWithoutThetas(uww_lr, th1, th2);
        vxy_lr = GetWithoutThetas(vww_lr, th1, th2);
        
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        if SETTINGS.PLOT_GRAPHS_LRA == true
        PlotSingularValues(fxy, gxy, fxy_lr, gxy_lr, fww, a_gww, fww_lr, a_gww_lr)
        end
        
        
        
    case 'None'
        % Dont Apply SNTLN improvements
        
        % Get preprocessed form
        fww = GetWithThetas(fxy, th1, th2);
        gww = GetWithThetas(gxy, th1, th2);
        a_gww = alpha.*gww;
        
        [uww, vww] = GetCofactors_Bivariate_2Polys(fww, a_gww, k1, k2);
        
        uxy_lr = GetWithoutThetas(uww, th1, th2);
        vxy_lr = GetWithoutThetas(vww, th1, th2);
        fxy_lr = fxy;
        gxy_lr = gxy;
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        % Plot Singular Values
        fww_lr = GetWithThetas(fxy_lr, th1_lr, th2_lr);
        a_gww_lr = alpha_lr.* GetWithThetas(gxy_lr, th1_lr, th2_lr);
        
        if SETTINGS.PLOT_GRAPHS_LRA == true
            PlotSingularValues(fxy, gxy, fxy_lr, gxy_lr, fww, a_gww, fww_lr, a_gww_lr);
        end
        
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
        
    otherwise
        error('bool_SNTLN is either (Standard SNTLN) or (None)')
end

end


function [] = PlotSingularValues(fxy, gxy, fxy_lr, gxy_lr, fww, a_gww, fww_lr, a_gww_lr)
%
% % Inputs
%
% fxy : (matrix)
% 
% gxy : (Matrix)
%
% fxy_lr : (Matrix)
%
% gxy_lr : (Matrix)
%
% fww : (Matrix)
%
% a_gww : (Matrix)
%
% fww_lr : (Matrix)
%
% gww_lr : (Matrix)


global SETTINGS


S1 = BuildDTQ_Bivariate_2Polys(fxy, gxy, 1, 1);
S2 = BuildDTQ_Bivariate_2Polys(fxy_lr, gxy_lr, 1, 1);
S3 = BuildDTQ_Bivariate_2Polys(fww, a_gww, 1, 1);
S4 = BuildDTQ_Bivariate_2Polys(fww_lr, a_gww_lr, 1, 1);

vSingularValues1 = svd(S1);
vSingularValues2 = svd(S2);
vSingularValues3 = svd(S3);
vSingularValues4 = svd(S4);

if (SETTINGS.PLOT_GRAPHS)
    figure()
    hold on
       
    plot(log10(vSingularValues1),'-s','DisplayName','$f(x,y), g(x,y)$');
    plot(log10(vSingularValues2),'-s','DisplayName','$f(x,y), g(x,y) LR$');
    
    plot(log10(vSingularValues3),'-s','DisplayName','$\tilde{f}(\omega_{1}, \omega_{2}), \alpha \tilde{g}(\omega_{1}, \omega_{2})$');
    
    plot(log10(vSingularValues4),'-s','DisplayName','$\tilde{f}(\omega_{1}, \omega_{2}), \alpha \tilde{g}(\omega_{1}, \omega_{2}) LR$');
    legend(gca,'show');
    hold off
end

% Legend and labels
l = legend(gca, 'show');
set(l, {'Interpreter'}, {'latex'})
set(l, {'FontSize'}, {20})
set(l, {'Location'}, {'southwest'})

xlabel('$i$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$\log_{10} \left( \sigma_{i} \right)$', 'Interpreter', 'latex', 'FontSize', 20)

grid on
box on



% Set location of window and size
m_left = 100;
m_bottom = 100;
m_width = 600;
m_height = 600;

set(gcf, 'Position', [m_left, m_bottom, m_width, m_height]);

% Position of figure within window
myplot = gca;
myval_side = 0.12;
myval_base = 0.10;
set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])


hold off

end

