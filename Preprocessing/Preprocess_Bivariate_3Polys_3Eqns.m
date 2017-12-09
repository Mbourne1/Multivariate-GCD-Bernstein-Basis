function [GM_fxy, GM_gxy, GM_hxy, alpha, beta, gamma, th1, th2] = ...
    Preprocess_Bivariate_3Polys_3Eqns(fxy, gxy, hxy, k1, k2)
%
% % Inputs
%
% fxy : (Matrix)
%
% gxy : (Matrix)
%
% hxy : (Matrix)
%
% k1 : (Int)
%
% k2 : (Int)
%
% % Outputs
%
% GM_fxy1 : (Float)
%
% GM_fxy2 : (Float)
%
% GM_gxy : (Float)
%
% GM_hxy : (Float)
%
% alpha : (Float)
%
% beta : (Float)
%
% gamma : (Float)
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}


% Get Degree of f(x,y), g(x,y) and h(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);


% Get Geometric Mean
GM_fxy = GetMean_2Partitions(fxy, n1, n2, o1, o2, k1, k2);
GM_gxy = GetMean_2Partitions(gxy, m1, m2, o1, o2, k1, k2);
GM_hxy = GetMean_2Partitions(hxy, m1, m2, n1, n2, k1, k2);


% Normalise by Geometric Mean
fxy_n = fxy ./ GM_fxy;
gxy_n = gxy ./ GM_gxy;
hxy_n = hxy ./ GM_hxy;



global SETTINGS

if (SETTINGS.BOOL_ALPHA_THETA)
    
    
    % %
    % Get optimal values of alpha and theta
    [max_matrix_fxy1, min_matrix_fxy1] = GetMaxMin(fxy_n, n1 - k1, n2 - k2);
    [max_matrix_fxy2, min_matrix_fxy2] = GetMaxMin(fxy_n, o1 - k1, o2 - k2);
    [max_matrix_gxy1, min_matrix_gxy1] = GetMaxMin(gxy_n, m1 - k1, m2 - k2);
    [max_matrix_gxy2, min_matrix_gxy2] = GetMaxMin(gxy_n, o1 - k1, o2 - k2);
    [max_matrix_hxy1, min_matrix_hxy1] = GetMaxMin(hxy_n, m1 - k1, m2 - k2);
    [max_matrix_hxy2, min_matrix_hxy2] = GetMaxMin(hxy_n, n1 - k1, n2 - k2);
    
    %alpha = OptimalAlpha(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g);
    
    [alpha, beta, gamma, th1, th2] = ...
        OptimalAlphaBetaGammaTheta_3Polys_3Eqns(...
        max_matrix_fxy1, min_matrix_fxy1, ...
        max_matrix_fxy2, min_matrix_fxy2, ...
        max_matrix_gxy1, min_matrix_gxy1,...
        max_matrix_gxy2, min_matrix_gxy2,...
        max_matrix_hxy1, min_matrix_hxy1,...
        max_matrix_hxy2, min_matrix_hxy2);
    
    
    
    
    
    
else
    
    alpha =1;
    beta = 1;
    gamma = 1;
    th1 = 1;
    th2 = 1;
    
end


end




function [] = PlotCoefficients(arrPolys, arrNames)
%
% % Inputs
%
% arrPolys : (Array of Matrices)
%
% arrNames : (Array of Strings)


nPolys = length(arrPolys);

arrColor = {'r','b'};

figure()
hold on

for i = 1:1:nPolys
    
    fxy = arrPolys{i};
    name = arrNames{i};
    
    [m1, m2] = GetDegree_Bivariate(fxy);
    
    nCoefficients = (m1 + 1) * (m2 + 1);
    
    v_fxy = GetAsVector_Version1((fxy), m1, m2);
    x_vec = 1:1:nCoefficients;
    
    plot(x_vec, log10(v_fxy), 'Color', arrColor{i}, 'DisplayName',name)
    
    % Get coefficient of maximum and minimum magnitude in f(x,y) and f(w,w)
    max_mag_fxy = max(log10(abs(v_fxy)));
    min_mag_fxy = min(log10(abs(v_fxy)));
    
    hline([max_mag_fxy, min_mag_fxy], {arrColor{i},arrColor{i}});
    
    xlim([1,nCoefficients]);
    
end


legend(gca,'show');
hold off


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


function [f_max,f_min,g_max,g_min] = GetMaxMinPairs(fxy, gxy, k1, k2)
% Get the maximum and minimum entries of each coefficient of f(x,y) and
% g(x,y) as appearing in the Sylvester subresultant S_{k_{1},k_{2}}(f,g)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of f(x,y)
%
% gxy : (Matrix) Coefficients of g(x,y)
%
% k1 : (Int)
%
% k2 : (Int)

% Get degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Get the maximum and minimum of each entry of f(x,y) in S_{k_{1},k_{2}}
[max_mtrx_f, min_mtrx_f] = GetMaxMin(fxy, n1 - k1, n2 - k2);

% Get the maximum and minimum of each entry of g(x,y) in S_{k_{1},k_{2}}
[max_mtrx_g, min_mtrx_g] = GetMaxMin(gxy, m1 - k1, m2 - k2);


f_max = max(max(max_mtrx_f));
f_min = min(min(min_mtrx_f));
g_max = max(max(max_mtrx_g));
g_min = min(min(min_mtrx_g));

end


