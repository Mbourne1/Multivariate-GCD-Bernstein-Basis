function [fxy_lr, gxy_lr, uxy_lr, vxy_lr, dxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    APF_Bivariate_2Polys(fxy, gxy, uxy, vxy, t1, t2, alpha, th1, th2)
% 
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y) 
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% uxy : (Matrix) Coefficients of the polynomial u(x,y) 
%
% vxy : (Matrix) Coefficients of the polynomial v(x,y)
%
% [t1, t2] : [(Int) (Int)] : Relative degrees of d(x,y)
%
% alpha : (Float): Optimal value of \alpha
%
% th1 th2 : [Float Float] : Optimal values of \theta_{1} and \theta_{2}
%
% % Outputs
%
% fxy_lr : (Matrix) Coefficients of f(x,y) after APF
%
% gxy_lr : (Matrix) Coefficients of g(x,y) after APF
%
% uxy_lr : (Matrix) Coefficients of u(x,y) after APF
%
% vxy_lr : (Matrix) Coefficients of v(x,y) after APF
%
% dxy_lr : (Matrix) Computed coefficients of d(x,y)
%
% alpha_lr : (Float) Optimal value of \alpha after APF
% 
% [th1_lr, th2_lr] : (Float Float) Optimal values of \theta_{1} and \theta_{2} after APF

global SETTINGS
switch SETTINGS.APF_METHOD
    
    
    case 'Standard APF Nonlinear'
        
        error([mfilename ' : This code branch is not complete']);
        
    case 'Standard APF Linear'
        
        error([mfilename ' : This code branch is not complete']);
        
    case 'None'
        
        % Get f(\omega_{1},\omega_{2}) and \alpha g(\omega_{1},\omega_{2})
        fww = GetWithThetas(fxy, th1, th2);
        a_gww = alpha.*GetWithThetas(gxy, th1, th2);
        
        % Get u(\omega_{1},\omega_{2}) and v(\omega_{1},\omega_{2})
        uww = GetWithThetas(uxy, th1, th2);
        vww = GetWithThetas(vxy, th1, th2);
        
        % Get d(\omega_{1},\omega_{2})
        [dww] = GetGCD_Coefficients_Bivariate_2Polys(uww, vww, fww, a_gww, t1, t2);
        
        % Get d(x,y)
        dxy_lr = GetWithoutThetas(dww, th1, th2);
        
        % Output are unchanged from input
        fxy_lr = fxy;
        gxy_lr = gxy;
        uxy_lr = uxy;
        vxy_lr = vxy;
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
        SETTINGS.APF_REQ_ITE = 0;
        
    otherwise
        
        error([ mfilename ' : Not valid APF Method']);
end


end