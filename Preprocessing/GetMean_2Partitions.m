function mean = GetMean_2Partitions(fxy, n1, n2, o1, o2, k1, k2)
% GetMean_2Partitions
%
% Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% n1 : (Int) 
% 
% n2 : (Int)
%
% o1 : (Int)
%
% o2 : (Int)
%
% k1 : (Int)
%
% k2 : (Int)



global SETTINGS

switch SETTINGS.MEAN_METHOD
    
    case 'Geometric Mean Matlab Method'
        mean = GetGeometricMean_2Partitions(fxy, n1, n2, o1, o2, k1, k2);
        
        
    case 'Arithmetic Mean'
        error('not completed')

    case 'None'
        mean = 1;
        
    otherwise
        error('error : Not a valid mean computation method')
end