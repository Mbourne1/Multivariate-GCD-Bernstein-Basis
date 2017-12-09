function [GM_fxy, GM_gxy, GM_hxy, lambda, mu, rho, th1, th2] = ...
    Preprocess_Bivariate_3Polys(fxy, gxy, hxy, k1, k2)


global SETTINGS

switch SETTINGS.SYLVESTER_MATRIX_3POLY_N_EQUATIONS
    
    
    case '2'
        [GM_fxy, GM_gxy, GM_hxy, lambda, mu, rho, th1, th2] = ...
            Preprocess_Bivariate_3Polys_2Eqns(fxy, gxy, hxy, k1, k2);
        
    case '3'
        
        [GM_fxy, GM_gxy, GM_hxy, lambda, mu, rho, th1, th2] = ...
            Preprocess_Bivariate_3Polys_3Eqns(fxy, gxy, hxy, k1, k2);
        
        
    otherwise
        error('err variable nEquations must be 2 or 3')
        
        
end



end