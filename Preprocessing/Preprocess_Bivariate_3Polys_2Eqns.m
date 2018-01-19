function [GM_fxy, GM_gxy, GM_hxy, alpha, beta, gamma, th1, th2] = ...
    Preprocess_Bivariate_3Polys_2Eqns(fxy, gxy, hxy, k1, k2)
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
GM_gxy = GetMean(gxy, m1 - k1, m2 - k2);
GM_hxy = GetMean(hxy, o1 - k1, o2 - k2);
    
    
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

    [max_matrix_gxy, min_matrix_gxy] = GetMaxMin(gxy_n, m1 - k1, m2 - k2);
    [max_matrix_hxy, min_matrix_hxy] = GetMaxMin(hxy_n, m1 - k1, m2 - k2);
    
    %alpha = OptimalAlpha(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g);
    
    [alpha, beta, gamma, th1, th2] = ...
        OptimalAlphaBetaGammaTheta_3Polys_2Eqns(...
        max_matrix_fxy1, min_matrix_fxy1, ...
        max_matrix_fxy2, min_matrix_fxy2, ...
        max_matrix_gxy, min_matrix_gxy,...
        max_matrix_hxy, min_matrix_hxy);
    

    
    
else
    
    alpha =1;
    beta = 1;
    gamma = 1;
    th1 = 1;
    th2 = 1;
    
end
    

end




function [] = PlotCoefficients(arrPolys,arrNames)
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




