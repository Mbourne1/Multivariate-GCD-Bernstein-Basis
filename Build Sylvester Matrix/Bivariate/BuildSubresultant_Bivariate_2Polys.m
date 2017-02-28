function [Sk1k2] = BuildSubresultant_Bivariate_2Polys(fxy, gxy, k1, k2)
%
% % Inputs
%
% fxy : Coefficients of the polynomial f(x,y)
%
% gxy : Coefficients of the polynomial g(x,y)
%
% k1 : 
%
% k2


% Get degree of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

global SETTINGS

switch SETTINGS.SYLVESTER_BUILD_METHOD
    case 'T'
        
        T1 = BuildT1_Bivariate(fxy, n1-k1, n2-k2);
        T2 = BuildT1_Bivariate(gxy, m1-k1, m2-k2);
        
        Sk1k2 = [T1 T2];
        
    case 'DT'
        
        BuildD_Bivariate_2Polys(m1, m2, n1-k1, n2-k2);
        T1 = BuildT1_Bivariate(fxy, n1-k1, n2-k2);
        T2 = BuildT1_Bivariate(gxy, m1-k1, m2-k2);
        
        Sk1k2 = D*[T1 T2];
        
    case 'DTQ'
        
        % Sk1k2 = BuildDTQ_Bivariate_2Polys(fxy, gxy, k1, k2);
        
        Sk1k2 = BuildDTQ_Bivariate_2Polys_NewMethod(fxy, gxy, k1, k2);
        
        
        
    case 'TQ'
        
        T1 = BuildT1_Bivariate(fxy, n1-k1, n2-k2);
        T2 = BuildT1_Bivariate(gxy, m1-k1, m2-k2);
        
        Q = BuildQ_Bivariate_2Polys(n1-k1,n2-k2,m1-k1,m2-k2);
        
        Sk1k2 = [T1 T2] * Q;
end