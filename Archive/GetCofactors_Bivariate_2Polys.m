function [uxy,vxy] = GetCofactors_Bivariate_2Polys(fxy, gxy, m, n, t, t1, t2, alpha, th1, th2)
%
% % Inputs
%
% [fxy, gxy] : Coefficients of the polynomials f(x,y) and g(x,y0
%
% [m, n] :
%
% t :
%
% [t1, t2] :
%
% alpha
%
% [th1, th2]
%
% % Outputs
%
% [uxy, vxy]



global SETTINGS

% Get f(\omega_{1},\omega_{2}) from f(x,y)
fww = GetWithThetas(fxy, th1, th2);

% Get \alpha * g(\omega_{1},\omega_{2}) from g(x,y)
gww = GetWithThetas(gxy, th1, th2);
a_gww = alpha.* gww;


switch SETTINGS.CALC_METHOD
    case 'Relative'
        
        [uww, vww] = GetCofactors_Relative_Bivariate_2Polys(fww, a_gww,t1,t2);
        
    case 'Total'
        
        [uww, vww] = GetCofactors_Total_Bivariate_2Polys(fww, a_gww, m, n, t);
        
    otherwise
        
        error([mfilename ':' 'Error degree Calc method is either (Relative) or (Total)']);
        
end


uxy = GetWithoutThetas(uww, th1, th2);
vxy = GetWithoutThetas(vww, th1, th2);


end