function [GM_fx, GM_gx, alpha, th1, th2] = Preprocess_Bivariate_2Polys(fxy, gxy, k1, k2)
% Preprocess(fxy,gxy,k1,k2)
%
% Get lamda,mu, alpha, theta_{1} and theta_{2}
%
% Inputs.
%
% [fxy, gxy] : Coefficients of f(x,y) and g(x,y)
%
% k1 : Degree of GCD d(x,y) with respect to x
%
% k2 : Degree of GCD d(x,y) with respect to y
%
% Outputs.
%
% [GM_fx, GM_gx] : Geometric Mean of entries of f(x,y) and g(x,y) in k-th subresultant matrix.
%
% alpha : Optimal value of alpha
%
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}

% Global variables
global SETTINGS

% Get degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Preproecessor One - Normalise by geometric mean
GM_fx = GetMean(fxy, n1-k1, n2-k2);

GM_gx = GetMean(gxy, m1-k1, m2-k2);


% Normalise f(x,y) by geometric mean
fxy_n = fxy./ GM_fx;

% Normalise g(x,y) by geometric mean
gxy_n = gxy./ GM_gx;

if (SETTINGS.BOOL_ALPHA_THETA)
    
    % %
    % Get maximum and minimum entries of f and g in the unprocessed
    % form without normalization by geometric mean.
    %[max_fxy, min_fxy, max_gxy, min_gxy] = GetMaxMinPairs(fxy,gxy,k1,k2);
    %PrintToFile(m1,m2,n1,n2,k1,k2,max_fxy,min_fxy,max_gxy,min_gxy,1,1,1,1,1);
    
    % %
    % Get Maximum and minimum Entries of f and g in the normalised but
    % unprocessed form.
    [max_fxy_n, min_fxy_n, max_gxy_n, min_gxy_n] = GetMaxMinPairs(fxy_n, gxy_n, k1, k2);
    
    PrintToFile(m1, m2, n1, n2, k1, k2, max_fxy_n, min_fxy_n, max_gxy_n, min_gxy_n, 1, 1, 1, GM_fx, GM_gx);
    
    
    % %
    % Get optimal values of alpha and theta
    [max_matrix_fxy, min_matrix_f] = GetMaxMin(fxy_n, n1-k1, n2-k2);
    [max_matrix_gxy, min_matrix_g] = GetMaxMin(gxy_n, m1-k1, m2-k2);
    
    %alpha = OptimalAlpha(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g);
    
    [alpha, th1, th2] = OptimalAlphaTheta(max_matrix_fxy, min_matrix_f, max_matrix_gxy, min_matrix_g);
    
    %fprintf('Optimal Value of alpha : %f \n',alpha)
    %fprintf('Optimal Value of theta1 : %f \n',th1)
    %fprintf('Optimal Value of theta2 : %f \n',th2)
    
    fww = GetWithThetas(fxy_n, th1, th2);
    gww = GetWithThetas(gxy_n, th1, th2);
    
    % %
    % Get Maximum and minimum entries of f and g in the normalised and
    % preprocessed form
    a_gww = alpha.*gww;
    
    if( SETTINGS.PLOT_GRAPHS)
        
        
        %                 figure('name','Preproc')
        %
        %                 subplot(1,2,1)
        %                 hold on
        %                 nCoefficients_fxy = (m1+1)*(m2+1);
        %                 nCoefficients_gxy = (n1+1)*(n2+1);
        %
        %                 xvals_fxy = 1:1:nCoefficients_fxy;
        %                 xvals_gxy = 1:1:nCoefficients_gxy;
        %
        %                 plot(xvals_fxy,log10(GetAsVector(fxy)),'DisplayName','f(x,y)');
        %                 plot(xvals_fxy,log10(GetAsVector(fww)),'DisplayName','f(\omega,\omega)');
        %                 legend(gca,'show');
        %                 hold off
        %
        %                 subplot(1,2,2)
        %                 hold on
        %                 plot(xvals_gxy,log10(GetAsVector(gxy)),'DisplayName','g(x,y)');
        %                 plot(xvals_gxy,log10(GetAsVector(a_gww)),'DisplayName','\alpha g(\omega,\omega)');
        %                 legend(gca,'show');
        %                 hold off
        
    end
    
    
    
    [max_fww,min_fww,max_gww,min_gww] = GetMaxMinPairs(fww, a_gww, k1, k2);
    
    PrintToFile(m1, m2, n1, n2, k1, k2, max_fww, min_fww, max_gww, min_gww, alpha, th1, th2, GM_fx, GM_gx);
    
    %fprintf([mfilename ' : ' sprintf('Condition S(f(x,y),g(x,y)) : %2.4f \n',cond(BuildDTQ(fxy,gxy,k1,k2)))]);
    %fprintf([mfilename ' : ' sprintf('Condition S(f(w),alpha.*g(w)) : %2.4f \n',cond(BuildDTQ(fww,alpha.*gww,k1,k2)))]);
    
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