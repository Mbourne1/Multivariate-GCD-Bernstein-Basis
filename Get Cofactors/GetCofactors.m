function [uxy,vxy] = GetCofactors(fxy,gxy,m,n,t,t1,t2,alpha,th1,th2)


global SETTINGS

% Get f(\omega_{1},\omega_{2}) from f(x,y)
fww = GetWithThetas(fxy,th1,th2);

% Get \alpha * g(\omega_{1},\omega_{2}) from g(x,y)
gww = GetWithThetas(gxy,th1,th2);
a_gww = alpha.* gww;


switch SETTINGS.CALC_METHOD
    case 'Relative'
        
        [uww, vww] = GetCofactors_Relative(fww, a_gww,t1,t2);
        
    case 'Total'
        
        [uww, vww] = GetCofactors_Total(fwww, alpha.*gww_n,m,n,t);
        
    otherwise
        
        error([mfilename ':' 'Error degree Calc method is either (Relative) or (Total)']);
        
end


uxy = GetWithoutThetas(uww,th1,th2);
vxy = GetWithoutThetas(vww,th1,th2);


end