function mean = GetMean(fxy,n1_k1,n2_k2)
% Get the mean of the coefficients of f(x,y) in the matrix
% C_{n1-k2,n2-k2}(f).

global SETTINGS

switch SETTINGS.MEAN_METHOD
    case 'Geometric Mean Matlab Method'
        mean = GetGeometricMean(fxy,n1_k1,n2_k2);
    otherwise
        error('err')
end