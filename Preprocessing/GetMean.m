function mean = GetMean(fxy, n1_k1, n2_k2)
% GetMean(fxy,n1_k1,n2_k2)
%
% Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% n1_k1 : (Int) Degree of polynomial v(x,y) with respect to x
%
% n2_k2 : (Int) Degree of polynomial v(x,y) with respect to y


% Get the mean of the coefficients of f(x,y) in the matrix
% C_{n1-k2,n2-k2}(f).

global SETTINGS

switch SETTINGS.MEAN_METHOD
    
    case 'Geometric Mean Matlab Method'
        mean = GetGeometricMean(fxy, n1_k1, n2_k2);
        
    case 'Geometric Mean My Method'
        mean = GetGeometricMean_MyMethod(fxy, n1_k1, n2_k2);
        
    case 'Arithmetic Mean'
        mean = GetArithmeticMean(fxy, n1_k1, n2_k2);

    case 'None'
        mean = 1;
        
    otherwise
        error('error : Not a valid mean computation method')
end