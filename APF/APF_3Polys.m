function [fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, dxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    APF_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, t, t1, t2, alpha, th1, th2)
%
% % Inputs
%
% [fxy, gxy, hxy] : Coefficients of the polynomials f(x,y), g(x,y) and
% h(x,y)
%
% [uxy, vxy, wxy] : Coefficients of the polynomials u(x,y), v(x,y) and
% w(x,y)
%
% [m, n, o] : Total degree of f(x,y), g(x,y) and h(x,y)
%
% [t , t1, t2] : Total and relative degrees of d(x,y)
%
% alpha : Optimal value of \alpha
%
% [th1 , th2] : Optimal values of \theta_{1} and \theta_{2}
%
% % Outputs
%
% [fxy_lr, gxy_lr, hxy_lr] : Coefficients of f(x,y), g(x,y) and h(x,y) after APF
%
% [uxy_lr, vxy_lr, wxy_lr] : Coefficients of u(x,y), v(x,y) and w(x,y) after APF
%
% dxy_lr : Computed coefficients of d(x,y)
%
% alpha_lr : Optimal value of \alpha after APF
% 
% [th1_lr, th2_lr] : Optimal values of \theta_{1} and \theta_{2} after APF

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
        hww = GetWithThetas(hxy, th1, th2);
        
        
        % Get u(\omega_{1},\omega_{2}) and v(\omega_{1},\omega_{2})
        uww = GetWithThetas(uxy, th1, th2);
        vww = GetWithThetas(vxy, th1, th2);
        www = GetWithThetas(wxy, th1, th2);
        
        % Get d(\omega_{1},\omega_{2})
        [dww] = GetGCD_Coefficients_3Polys(fww, a_gww, hww, uww, vww,www, m, n, o, t, t1, t2);
        
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