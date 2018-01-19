function [Sk1k2] = BuildSubresultant_Bivariate_3Polys_3Eqns(fxy, gxy, hxy, k1, k2)
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
% k1 : (Int) Index of subresultant matrix
%
% k2 : (Int) Index of subresultant matrix


% Get degree of polynomial f(x,y), g(x,y) and h(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);

global SETTINGS

% TO DO - Write code for 'T Version 2' 'DT Version 2' 'TQ Version 2'


switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    
    case 'T'
        
        T = BuildT_Bivariate_3Polys_3Eqns(fxy, gxy, hxy, k1, k2);
        
        Sk1k2 = T;
        
    case 'DT'
        
        D = BuildD_Bivariate_3Polys_3Eqns(m1, m2, n1, n2, o1, o2, k1, k2);
        T = BuildT_Bivariate_3Polys_3Eqns(fxy, gxy, hxy, k1, k2);
        
        Sk1k2 = D * T;
        
        
    case 'TQ'

        T = BuildT_Bivariate_3Polys_3Eqns(fxy, gxy, hxy, k1, k2);
        Q = BuildQ_Bivariate_3Polys(m1, m2, n1, n2, o1, o2, k1, k2);

        Sk1k2 = T * Q;
        
    case 'DTQ'
        
        D = BuildD_Bivariate_3Polys_3Eqns(m1, m2, n1, n2, o1, o2, k1, k2);
        T = BuildT_Bivariate_3Polys_3Eqns(fxy, gxy, hxy, k1, k2);
        Q = BuildQ_Bivariate_3Polys(m1, m2, n1, n2, o1, o2, k1, k2);

        
        Sk1k2 = D * T * Q;      
        
        
        
    otherwise
        error('Not a valid Sylvester_Build Method');
    
end