function [t1, t2, GM_fxy, GM_gxy, GM_hxy, lambda, mu, rho, th1, th2] = GetGCDDegree_Bivariate_3Polys_Linear(fxy, gxy, hxy, t_known, str_variable)
%
%
% % Inputs
%
% fxy : (Matrix) Coefficients of f(x,y)
%
% gxy : (Matrix) Coefficients of g(x,y)
%
% t_known : (Int) Either degree of GCD wrt x or wrt y
%
% str_variable : (String) 'x' or 'y' denotes which variable t is known for.
%
%
% % Outputs
%
% t1 : (Int) The degree of the GCD with respect to x
%
% t2 : (Int) The degree of the GCD with respect to y
%
% GM_fxy : (Float) Geometric mean of the coefficients of f(x,y)
%
% GM_gxy : (Float) Geometric mean of the coefficients of g(x,y)
%
% GM_hxy : (Float) Geometric mean of the coefficients of h(x,y)
%
% lambda : (Float) Optimal value \lambda
%
% mu : (Float) Optimal value \mu
%
% rho : (Float) Optimal value \rho
%
% th1 : (Float) Optimal value \theta_{1}
%
% th2 : (Float) Optimal value \theta_{2}

% Get the degree of the polynomials f(x,y), g(x,y) and h(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);


if strcmp(str_variable , 'x')
    
    k1 = t_known;
    t1 = t_known;
    
    % t2 is to be determined
    
    % Get limits for k2
    lowerLimit_k = 1;
    upperLimit_k = min([m2, n2, o2]);
    limits_k = [lowerLimit_k upperLimit_k];
    
    
    nSubresultants = upperLimit_k - lowerLimit_k + 1;
    
else
    k2 = t_known;
    t2 = t_known;
    % t1 is to be determined
    
    
    lowerLimit_k = 1;
    upperLimit_k = min([m1, n1, o1]);
    limits_k = [lowerLimit_k, upperLimit_k];
    
    % Get number of subresultant matrices
    nSubresultants = upperLimit_k - lowerLimit_k + 1;
    
end


% Initialise matrices to store values
vLambda = zeros(nSubresultants);
vMu = zeros(nSubresultants);
vRho = zeros(nSubresultants);

vTheta1 = zeros(nSubresultants);
vTheta2 = zeros(nSubresultants);

vGM_fxy = zeros(nSubresultants);
vGM_gxy = zeros(nSubresultants);
vGM_hxy = zeros(nSubresultants);

arr_Sk1k2 = cell(nSubresultants);
arr_SingularValues = cell(nSubresultants);



for i = 1 : 1 : nSubresultants
    
    
    if strcmp(str_variable, 'x')
        
        k2 = lowerLimit_k + (i - 1);
        
    else
        
        k1 = lowerLimit_k + (i - 1);
        
    end
    
    
    
    
    % Preprocessing
    [GM_fxy, GM_gxy, GM_hxy, lambda, mu, rho, th1, th2] = ...
        Preprocess_Bivariate_3Polys(fxy, gxy, hxy, k1, k2);
    
    % Divide f(x,y), g(x,y) and h(x,y) by their respective geometric mean
    fxy_matrix_n = fxy ./ GM_fxy;
    gxy_matrix_n = gxy ./ GM_gxy;
    hxy_matrix_n = hxy ./ GM_hxy;
    
    
    % Get f(w,w) from f(x,y) and g(w,w) from g(x,y)
    fww_matrix = GetWithThetas(fxy_matrix_n, th1, th2);
    gww_matrix = GetWithThetas(gxy_matrix_n, th1, th2);
    hww_matrix = GetWithThetas(hxy_matrix_n, th1, th2);
    
    
    % Build the Sylvester subresultant matrix S_{k1,k2}(f,g)
    arr_Sk1k2{i,1} = BuildSubresultant_Bivariate_3Polys(lambda.*fww_matrix, mu.*gww_matrix, rho.*hww_matrix, k1, k2);
    
    arr_SingularValues{i,1} = svd(arr_Sk1k2{i});
    
    % Store the optimal alpha, th1, th2, gm_fx, gm_gx
    vLambda(i, 1) = lambda;
    vMu(i, 1) = lambda;
    vRho(i, 1) = rho;
    
    vTheta1(i, 1) = th1;
    vTheta2(i, 1) = th2;
    
    vGM_fxy(i, 1) = GM_fxy;
    vGM_gxy(i, 1) = GM_gxy;
    vGM_hxy(i, 1) = GM_hxy;
    
end




global SETTINGS


switch SETTINGS.RANK_REVEALING_METRIC
    
    case 'Minimum Singular Values'
        
        % Initialise vector to store minimum singular values
        vMinimumSingularValues = zeros(nSubresultants, 1);
        
        for i = 1:1:nSubresultants
            
            % Get SVD of unproc and processed Sylvester Surbesultant S_{k,k}
            vMinimumSingularValues(i) = min(arr_SingularValues{i});
            
        end
        
        % Set metric
        metric = vMinimumSingularValues;
        
        % Plot Graphs
        if(SETTINGS.PLOT_GRAPHS)
            %plotSingularValues_1Dimensional(arr_SingularValues, limits_k, limits_t);
            limits_t = limits_k;
            
            
            if strcmp(str_variable , 'x')
                % t2 is to be determined
                txtXlabel = "k_{2}";
                txtYlabel = "\sigma_{t_{1},k_{2}}";
                txtTitle = "Singular Values of St1k2";
                
            else
                % t1 is to be determined
                txtXlabel = "k_{1}";
                txtYlabel = "\sigma_{k_{1},t_{2}}";
                txtTitle = "Singular Values of Sk1t2";
                
                
            end
            
            plotMinimumSingularValues_1Dimensional(vMinimumSingularValues, ...
                limits_k, limits_t, ...
                txtTitle, txtXlabel, txtYlabel);
        end
        
    case 'R1 Row Norms'
        
        % Initialise vectors to store Max and Min Row Norms of R1_{k}
        arr_R1_RowNorm = cell(nSubresultants,1);
        vMaxRowNorm = zeros(nSubresultants,1);
        vMinRowNorm = zeros(nSubresultants,1);
        
        
        for i = 1 : 1 : nSubresultants
            
            % Get Norms of each row in the matrix R1
            arr_R1_RowNorm{i} = sqrt(sum(arr_R1{i}.^2,2))./norm(arr_R1{i});
            vMaxRowNorm = max(arr_R1_RowNorm{i});
            vMinRowNorm = min(arr_R1_RowNorm{i});
            
        end
        
        % Set metric
        metric = vMinRowNorm./ vMaxRowNorm;
        
        % Plot graphs
        
    case 'R1 Row Diagonals'
        
        % Initialise vectors
        arr_R1_diag = cell(nSubresultants,1);
        vMaxDiagonal_R1 = zeros(nSubresultants,1);
        vMinDiagonal_R1 = zeros(nSubresultants,1);
        
        % Get ONLY the diagonal elements and normalise them.
        for i = 1:1:nSubresultants
            arr_R1_diag{i} = diag(arr_R1{i});
            vMaxDiagonal_R1(i) = max(abs(arr_R1_diag{i}));
            vMinDiagonal_R1(i) = min(abs(arr_R1_diag{i}));
        end
        
        % Set metric
        metric = vMinDiagonal_R1 ./ vMaxDiagonal_R1;
        
        % Plot graphs
        if(SETTINGS.PLOT_GRAPHS)
            %plotMaxMinDiagonals_R1_Total(vMaxDiagonal_R1, vMinDiagonal_R1, myLimits, limits_t);
            %plotDiagonals_R1(arr_R1_diag);
        end
        
    case 'Residuals'
        error('Code not complete')
    otherwise
        error('err');
        
end



if upperLimit_k == lowerLimit_k
    
    t_star = GetGCDDegree_OneSubresultant(metric, limits_k);
    
    
    
else
    
    t_star = GetGCDDegree_MultipleSubresultants(metric, limits_k);
    
    
end

if t_star == 0
    th1 = 1;
    th2 = 1;
    lambda = 1;
    return;
end


if strcmp(str_variable, 'x')
    
    t2 = t_star;
    
else
    
    t1 = t_star;
    
end











end