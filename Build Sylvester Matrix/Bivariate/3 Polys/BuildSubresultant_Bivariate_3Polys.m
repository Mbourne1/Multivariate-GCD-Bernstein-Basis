function [Sk1k2] = BuildSubresultant_Bivariate_3Polys(fxy, gxy, hxy, k1, k2)
% Build the Sylvester subresultant matrix S_{k_{1},k_{2}}
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of the polynomial h(x,y)
%
% k1 : (Int) Index of subresultant matrix S_{k_{1}, k_{2}}(f, g, h)
%
% k2 : (Int) Index of subresultant matrix S_{k_{1}, k_{2}}(f, g, h)



global SETTINGS


switch SETTINGS.SYLVESTER_MATRIX_3POLY_N_EQUATIONS
    
    case '2'
        Sk1k2 = BuildSubresultant_Bivariate_3Polys_2Eqns(fxy, gxy, hxy, k1, k2);
        
    case '3'
        Sk1k2 = BuildSubresultant_Bivariate_3Polys_3Eqns(fxy, gxy, hxy, k1, k2);
        
    otherwise
        error('Not a valid Sylvester_Build Method');
    
end


end