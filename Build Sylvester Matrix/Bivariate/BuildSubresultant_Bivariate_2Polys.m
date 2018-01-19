function [Sk1k2] = BuildSubresultant_Bivariate_2Polys(fxy, gxy, k1, k2)
% Build the Sylvester subresultant matrix S_{k_{1},k_{2}}
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% k1 : (Int)
%
% k2 : (Int)


% Get degree of polynomial f(x,y) and g(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

global SETTINGS

% TO DO - Write code for 'T Version 2' 'DT Version 2' 'TQ Version 2'


switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    
    case 'T'
        
        T1 = BuildT1_Bivariate(fxy, n1 - k1, n2 - k2);
        T2 = BuildT1_Bivariate(gxy, m1 - k1, m2 - k2);
        
        Sk1k2 = [T1 T2];
        
    case 'DT'
        
        D = BuildD_Bivariate_2Polys(m1, m2, n1 - k1, n2 - k2);
        T1 = BuildT1_Bivariate(fxy, n1 - k1, n2 - k2);
        T2 = BuildT1_Bivariate(gxy, m1 - k1, m2 - k2);
        
        Sk1k2 = D*[T1 T2];
        
    case 'DTQ'
        
        D = BuildD_Bivariate_2Polys(m1, m2, n1 - k1, n2 - k2);
        T1 = BuildT1_Bivariate(fxy, n1 - k1, n2 - k2);
        T2 = BuildT1_Bivariate(gxy, m1 - k1, m2 - k2);
        Q = BuildQ_Bivariate_2Polys(n1 - k1, n2 - k2, m1 - k1, m2 - k2);
        
        Sk1k2 = D*[T1 T2]*Q;
        
    case 'DTQ Version 2'
        
        Sk1k2 = BuildDTQ_Bivariate_2Polys_Version2(fxy, gxy, k1, k2);
        
    case 'TQ'
        
        T1 = BuildT1_Bivariate(fxy, n1 - k1, n2 - k2);
        T2 = BuildT1_Bivariate(gxy, m1 - k1, m2 - k2);
        
        Q = BuildQ_Bivariate_2Polys(n1 - k1, n2 - k2, m1 - k1, m2 - k2);
        
        Sk1k2 = [T1 T2] * Q;
        
    case 'DTQ Denominator Removed'
        
        D = BuildD_Bivariate_2Polys(m1, m2, n1 - k1, n2 - k2);
        
        T1 = BuildT1_Bivariate(fxy, n1 - k1, n2 - k2) .* ...
            nchoosek(m1 + n1 - k1, n1 - k1) * nchoosek(m2 + n2 - k2, n2-k2);
        
        T2 = BuildT1_Bivariate(gxy, m1 - k1, m2 - k2) .* ...
            nchoosek(m1 + n1 - k1, m1 - k1) * nchoosek(m2 + n2 - k2, m2 - k2);
        
        Q = BuildQ_Bivariate_2Polys(n1 - k1, n2 - k2, m1 - k1, m2 - k2);
        
        Sk1k2 = D*[T1 T2]*Q;
        
        
    otherwise
        error('Not a valid Sylvester_Build Method');
    
end