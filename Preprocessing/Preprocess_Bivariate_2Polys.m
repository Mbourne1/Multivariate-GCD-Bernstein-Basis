function [GM_fx, GM_gx, alpha, th1, th2] = Preprocess_Bivariate_2Polys(fxy, gxy, k1, k2)
% Preprocess(fxy,gxy,k1,k2)
%
% Get lamda,mu, alpha, theta_{1} and theta_{2}
%
% Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% k1 : (Int) Degree of GCD d(x,y) with respect to x
%
% k2 : (Int) Degree of GCD d(x,y) with respect to y
%
% Outputs.
%
% GM_fx : (Float) Geometric Mean of entries of f(x,y) in k-th subresultant matrix.
%
% GM_gx : (Float) Geometric Mean of entries of g(x,y) in k-th subresultant matrix.
%
% alpha : (Float) Optimal value of alpha
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}

% Global variables
global SETTINGS

% Get degree of f(x,y) and g(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Preproecessor One - Normalise by geometric mean
GM_fx = GetMean(fxy, n1-k1, n2-k2);
GM_gx = GetMean(gxy, m1-k1, m2-k2);


% Normalise f(x,y) and g(x,y) by respective geometric means
fxy_n = fxy./ GM_fx;
gxy_n = gxy./ GM_gx;


if (SETTINGS.BOOL_ALPHA_THETA)
    
    
    % %
    % Get optimal values of alpha and theta
    [max_matrix_fxy, min_matrix_f] = GetMaxMin(fxy_n, n1-k1, n2-k2);
    [max_matrix_gxy, min_matrix_g] = GetMaxMin(gxy_n, m1-k1, m2-k2);
    
    %alpha = OptimalAlpha(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g);
    
    [alpha, th1, th2] = OptimalAlphaTheta(max_matrix_fxy, min_matrix_f, max_matrix_gxy, min_matrix_g);
    
    fww = GetWithThetas(fxy_n, th1, th2);
    gww = GetWithThetas(gxy_n, th1, th2);
    
    % %
    % Get Maximum and minimum entries of f and g in the normalised and
    % preprocessed form
    a_gww = alpha.*gww;
  
    [max_fww,min_fww,max_gww,min_gww] = GetMaxMinPairs(fww, a_gww, k1, k2);
    
    PrintToFile(m1, m2, n1, n2, k1, k2, max_fww, min_fww, max_gww, min_gww, alpha, th1, th2, GM_fx, GM_gx);

else
    
    alpha =1;
    th1 = 1;
    th2 = 1;
    
end
end

function [] = PrintToFile(m1,m2,n1,n2,k1,k2,max_f,min_f,max_g,min_g,alpha,th1,th2,lambda,mu)


global SETTINGS


fullFileName = 'Results/Results_Preprocessing.txt';


if exist('Results/Results_Preprocessing.txt', 'file')
    fileID = fopen('Results/Results_Preprocessing.txt','a');
    fprintf(...
        fileID,...
        '%s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s,  \t %s \n',...
        datetime(),...
        SETTINGS.EX_NUM,...
        num2str(m1),...
        num2str(m2),...
        num2str(n1),...
        num2str(n2),...
        num2str(k1),...
        num2str(k2),...
        num2str(max_f),...
        num2str(min_f),...
        num2str(max_g),...
        num2str(min_g),...
        num2str(alpha),...
        num2str(th1),...
        num2str(th2),...
        num2str(SETTINGS.BOOL_ALPHA_THETA),...
        SETTINGS.EMIN...
        );
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end


end


function [f_max,f_min,g_max,g_min] = GetMaxMinPairs(fxy,gxy,k1,k2)
% Get the maximum and minimum entries of each coefficient of f(x,y) and
% g(x,y) as appearing in the Sylvester subresultant S_{k_{1},k_{2}}(f,g)


% Get degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Get the maximum and minimum of each entry of f(x,y) in S_{k_{1},k_{2}}
[max_mtrx_f, min_mtrx_f] = GetMaxMin(fxy, n1-k1, n2-k2);

% Get the maximum and minimum of each entry of g(x,y) in S_{k_{1},k_{2}}
[max_mtrx_g, min_mtrx_g] = GetMaxMin(gxy, m1-k1, m2-k2);


f_max = max(max(max_mtrx_f));
f_min = min(min(min_mtrx_f));
g_max = max(max(max_mtrx_g));
g_min = min(min(min_mtrx_g));

end