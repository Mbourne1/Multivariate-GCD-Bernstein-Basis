function [lambda, mu, alpha, th1,th2] = Preprocess(fxy_matrix_delv,gxy_matrix_delv,k1,k2)
% Get lamda,mu, alpha, theta_{1} and theta_{2}
%
% Inputs.
%
% fxy_matrix_delv : Coefficients of f(x,y) in degree elevated form
%
% gxy_matrix_delv : Coefficients of g(x,y) in degree elevated form
%
% k : Degree of GCD d(x,y).
% 
% Outputs.
%
% lamda : Geometric Mean of entries of f(x,y) in k-th subresultant matrix.
%
% mu : Geometric Mean of entries of g(x,y) in k-th subresultant matrix.
%
% alpha : Optimal value of alpha
% 
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}

% Global variables
global BOOL_PREPROC

% Get degree of f(x,y)
[m1,m2] = GetDegree(fxy_matrix_delv);

% Get degree of g(x,y)
[n1,n2] = GetDegree(gxy_matrix_delv);


switch BOOL_PREPROC
    case 'y'
        
        % Preproecessor One - Normalise by geometric mean
        [lambda, mu] = GetGeometricMean(fxy_matrix_delv,gxy_matrix_delv,k1,k2);
        
        % Normalise f(x,y)
        fxy_matrix_n = fxy_matrix_delv./lambda;
        
        % Normalise g(x,y)
        gxy_matrix_n = gxy_matrix_delv./mu;
        
        % Get the maximum and minimum entries of f(x,y) in the
        % Sylvester matrix S(f,g)
        [max_mtrx_f, min_mtrx_f] = GetMaxMin(fxy_matrix_n,n1,n2,k1,k2);
        
        % Get the maximum and minimum entries of g(x,y) in the
        % Sylvester matrix S(f,g)
        [max_mtrx_g, min_mtrx_g] = GetMaxMin(gxy_matrix_n,m1,m2,k1,k2);
        
        % Get optimal values of alpha and theta
        [alpha, th1, th2] = OptimalAlphaTheta(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g);
        
    case 'n'
        alpha = 1;
        th1 = 1;
        th2 = 1;
        lambda = 1;
        mu = 1;
    otherwise
        error('err')
end