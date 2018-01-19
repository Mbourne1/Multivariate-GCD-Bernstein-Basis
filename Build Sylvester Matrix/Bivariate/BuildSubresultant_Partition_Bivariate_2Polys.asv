function [C1] = BuildSubresultant_Partition_Bivariate_2Polys(fxy, n1_k1, n2_k2)
% Build the Sylvester subresultant matrix S_{k_{1},k_{2}}
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% n1_k1 : (Int)
%
% n2_k2 : (Int)
%
% % Outputs
% 
% C1 : (Matrix) Subresutlant matrix partition

% Get degree of polynomial f(x,y) and g(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);

global SETTINGS

% TO DO - Write code for 'T Version 2' 'DT Version 2' 'TQ Version 2'


switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    
    case 'T'
        
        T1 = BuildT1_Bivariate(fxy, n1_k1, n2_k2);
        
        
        C1 = T1;
        
    case 'DT'
        
        D = BuildD_Bivariate_2Polys(m1, m2, n1_k1, n2_k2);
        T1 = BuildT1_Bivariate(fxy, n1_k1, n2_k2);
        
        C1 = D*T1;
        
    case 'DTQ'
        
        D = BuildD_Bivariate_2Polys(m1, m2, n1_k1, n2_k2);
        T1 = BuildT1_Bivariate(fxy, n1_k1, n2_k2);
        Q1 = BuildQ1_Bivariate(n1_k1, n2_k2);
        
        C1 = D*T1*Q1;
        
    case 'TQ'
        
        T1 = BuildT1_Bivariate(fxy, n1_k1, n2_k2);
        Q1 = BuildQ1_Bivariate(n1_k1, n2_k2);
        
        C1 = T1*Q1;
        
    case 'DTQ Denominator Removed'
        
        D = BuildD_Bivariate_2Polys(m1, m2, n1_k1, n2_k2);
        
        T1 = BuildT1_Bivariate(fxy, n1_k1, n2_k2) .* ...
            nchoosek(m1 + n1_k1, n1_k1) * nchoosek(m2 + n2_k2, n2_k2);
        
        Q1 = BuildQ1_Bivariate(n1_k1, n2_k2);
        
        C1 = D*T1*Q1;
        
        
    otherwise
        error('Not a valid Sylvester_Build Method');
    
end