function [uxy, vxy, wxy] = GetCofactors_3Polys(fxy, gxy, hxy, m, n, o, t, t1, t2, alpha, th1, th2)
%
% % Inputs
%
% [fxy, gxy, hxy] : Coefficients of the polynomials f(x,y), g(x,y) and
% h(x,y)
%
% [m, n, o] : Total degree of f(x,y) g(x,y) and h(x,y)
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
% [uxy, vxy, wxy] : Coefficients of u(x,y) v(x,y) and w(x,y)



global SETTINGS

% Get f(\omega_{1},\omega_{2}) from f(x,y)
fww = GetWithThetas(fxy, th1, th2);

% Get \alpha * g(\omega_{1},\omega_{2}) from g(x,y)
gww = GetWithThetas(gxy, th1, th2);
a_gww = alpha.* gww;

hww = GetWithThetas(hxy, th1, th2);


switch SETTINGS.CALC_METHOD
    case 'Relative'
        
        [uww, vww, www] = GetCofactors_Relative_3Polys(fww, a_gww, hww, t1, t2);
        
    case 'Total'
        
        [uww, vww, www] = GetCofactors_Total_3Polys(fww, alpha.*gww, hww, m, n, o, t);
        
    otherwise
        
        error([mfilename ':' 'Error degree Calc method is either (Relative) or (Total)']);
        
end


uxy = GetWithoutThetas(uww, th1, th2);
vxy = GetWithoutThetas(vww, th1, th2);
wxy = GetWithoutThetas(www, th1, th2);

end