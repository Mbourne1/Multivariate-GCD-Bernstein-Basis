function [lambda, mu, alpha, th1,th2] = Preprocess(fxy,gxy,k1,k2)
% Get lamda,mu, alpha, theta_{1} and theta_{2}
%
% Inputs.
%
% fxy_matrix_delv : Coefficients of f(x,y) in degree elevated form
%
% gxy_matrix_delv : Coefficients of g(x,y) in degree elevated form
%
% k1 : Degree of GCD d(x,y) with respect to x
%
% k2 : Degree of GCD d(x,y) with respect to y
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
global BOOL_ALPHA_THETA


% Get degree of f(x,y)
[m1,m2] = GetDegree(fxy);

% Get degree of g(x,y)
[n1,n2] = GetDegree(gxy);

% Preproecessor One - Normalise by geometric mean
lambda = GetMean(fxy,n1-k1,n2-k2);

mu = GetMean(gxy,n1-k1,n2-k2);

% Normalise f(x,y) by geometric mean
fxy_n = fxy./lambda;

% Normalise g(x,y) by geometric mean
gxy_n = gxy./mu;

switch BOOL_ALPHA_THETA
    case 'y'
        
        
        % %
        % Get maximum and minimum entries of f and g in the unprocessed
        % form without normalization by geometric mean.
        [max_fxy, min_fxy, max_gxy, min_gxy] = GetMaxMinPairs(fxy,gxy,k1,k2);
        %PrintToFile(m1,m2,n1,n2,k1,k2,max_fxy,min_fxy,max_gxy,min_gxy,1,1,1,1,1);

        % % 
        % Get Maximum and minimum Entries of f and g in the normalised but
        % unprocessed form.
        [max_fxy_n, min_fxy_n, max_gxy_n, min_gxy_n] = GetMaxMinPairs(fxy_n, gxy_n,k1,k2);
        %PrintToFile(m1,m2,n1,n2,k1,k2,max_fxy_n,min_fxy_n,max_gxy_n,min_gxy_n,1,1,1,lambda,mu);
        
        
        % %
        % Get optimal values of alpha and theta
        [max_mtrx_f, min_mtrx_f] = GetMaxMin(fxy_n,n1-k1,n2-k2);
        [max_mtrx_g, min_mtrx_g] = GetMaxMin(gxy_n,m1-k1,m2-k2);
        
        [alpha, th1, th2] = OptimalAlphaTheta(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g);
        
        fww = GetWithThetas(fxy_n,th1,th2);
        gww = GetWithThetas(gxy_n,th1,th2);
        
        % % 
        % Get Maximum and minimum entries of f and g in the normalised and
        % preprocessed form
        [max_fww,min_fww,max_gww,min_gww] = GetMaxMinPairs(fww,gww,k1,k2);
        %PrintToFile(m1,m2,n1,n2,k1,k2,max_fww,min_fww,max_gww,min_gww,alpha,th1,th2,lambda,mu);
        
    case 'n'
        alpha = 1;
        th1 = 1;
        th2 = 1;
    otherwise
        error('err')
end
end

function [] = PrintToFile(m1,m2,n1,n2,k1,k2,max_f,min_f,max_g,min_g,alpha,th1,th2,lambda,mu)


data = [m1 m2 n1 n2 k1 k2 max_f min_f max_g min_g alpha th1 th2 lambda mu];

dlmwrite('Preprocessing.txt',data,'-append',...
'delimiter','\t','roffset',0)



end


function [f_max,f_min,g_max,g_min] = GetMaxMinPairs(fxy,gxy,k1,k2)

% Get Degree of f(x,y)
[m1,m2] = GetDegree(fxy);
[n1,n2] = GetDegree(gxy);

[max_mtrx_f, min_mtrx_f] = GetMaxMin(fxy,n1-k1,n2-k2);
[max_mtrx_g, min_mtrx_g] = GetMaxMin(gxy,m1-k1,m2-k2);

f_max = max(max(max_mtrx_f));
f_min = min(min(min_mtrx_f));
g_max = max(max(max_mtrx_g));
g_min = min(min(min_mtrx_g));

end