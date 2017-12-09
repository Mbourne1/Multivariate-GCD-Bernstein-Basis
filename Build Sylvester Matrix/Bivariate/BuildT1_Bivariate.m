function T1 = BuildT1_Bivariate(fxy, n1_k1, n2_k2)
% BuildT1_Bivariate(fxy_matrix, n1_k1, n2_k2)
%
% Build the matrix T1, the first partition of the Sylvester matrix
% S_{k_{1},k_{2}}. T_{n_{1}-k_{1},n_{2}-k_{2}}(f) * [v] gives a vector of coefficients of
% the polynomial multiplication f(x,y) * v(x,y).
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% n1_k1 : (Int) Degree of polynomial v(x,y) with respect to x
%
% n2_k2 : (Int) Degree of polynomial v(x,y) with repsect to y
%
% % Outputs
%
% T1 : (Matrix) T_{n1-k1,n2-k2}(f(x,y))


global SETTINGS
switch SETTINGS.VECTORISATION_METHOD
    case 'Version 1'
        
        T1 = BuildT1_Bivariate_Version1(fxy, n1_k1, n2_k2);
        
    case 'Version 2'
        
        T1 = BuildT1_Bivariate_Version2(fxy, n1_k1, n2_k2);
        
    otherwise
        error('error')
        
end

end