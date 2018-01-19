function [t1, t2, GM_fx, GM_gx, alpha, th1, th2] = GetGCDDegree_Bivariate_2Polys_Linear(fxy, gxy, t_known, str_variable)
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
% t1 : (Int)
%
% t2 : (Int)
%
% GM_fx : (Float)
%
% GM_gx : (Float)
%
% alpha : (Float)
%
% th1 : (Float)
%
% th2 : (Float)


[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);



if strcmp(str_variable , 'x')
    
    k1 = t_known;
    t1 = t_known;
    % t2 is to be determined
    
    % Get limits for k2
    lowerLimit_k = 1;
    upperLimit_k = min(m2, n2);
    limits_k = [lowerLimit_k upperLimit_k];
    
    
    nSubresultants = upperLimit_k - lowerLimit_k + 1;
    
else
    k2 = t_known;
    t2 = t_known;
    % t1 is to be determined
    
    
    lowerLimit_k = 1;
    upperLimit_k = min(m1, n1);
    limits_k = [lowerLimit_k, upperLimit_k];
    
    nSubresultants = upperLimit_k - lowerLimit_k + 1;
    
end


% Initialise matrices to store values
vAlpha = zeros(nSubresultants);
vTheta1 = zeros(nSubresultants);
vTheta2 = zeros(nSubresultants);
vGM_fx = zeros(nSubresultants);
vGM_gx = zeros(nSubresultants);

arr_Sk1k2 = cell(nSubresultants);
arr_SingularValues = cell(nSubresultants);



for i = 1 : 1 : nSubresultants
    
    
    if strcmp(str_variable, 'x')
        
        k2 = lowerLimit_k + (i - 1);
        
    else
        
        k1 = lowerLimit_k + (i - 1);
        
    end
    
    
    % Preprocessing
    [GM_fx, GM_gx, alpha, th1, th2] = Preprocess_Bivariate_2Polys(fxy, gxy, k1, k2);
    
    % Divide f(x,y) and g(x,y) by geometric mean
    fxy_matrix_n = fxy ./ GM_fx;
    gxy_matrix_n = gxy ./ GM_gx;
    
    % Get f(w,w) from f(x,y) and g(w,w) from g(x,y)
    fww_matrix = GetWithThetas(fxy_matrix_n, th1, th2);
    gww_matrix = GetWithThetas(gxy_matrix_n, th1, th2);
    
    % Build the Sylvester subresultant matrix S_{k1,k2}(f,g)
    arr_Sk1k2{i,1} = BuildSubresultant_Bivariate_2Polys(fww_matrix, alpha.*gww_matrix, k1, k2);
    
    arr_SingularValues{i,1} = svd(arr_Sk1k2{i});
    
    % Store the optimal alpha, th1, th2, gm_fx, gm_gx
    vAlpha(i, 1) = alpha;
    vTheta1(i, 1) = th1;
    vTheta2(i, 1) = th2;
    vGM_fx(i, 1) = GM_fx;
    vGM_gx(i, 1) = GM_gx;
    

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
            
            limits_t = limits_k;
            plotSingularValues_1Dimensional(arr_SingularValues, limits_k, limits_t);
            plotMinimumSingularValues_1Dimensional(vMinimumSingularValues, limits_k, limits_t);
        end
        
    case 'Normalised Minimum Singular Values'
        
        % Initialise vector to store minimum singular values
        vMinimumSingularValues = zeros(nSubresultants, 1);
        
        for i = 1 : 1 : nSubresultants
            
            % Get SVD of unproc and processed Sylvester Surbesultant S_{k,k}
            vMinimumSingularValues(i) = min(arr_SingularValues{i});
            
        end
        
        % Set metric
        metric = vMinimumSingularValues ./ vMinimumSingularValues(1);
        
        % Plot Graphs
        if(SETTINGS.PLOT_GRAPHS)
            
            limits_t = limits_k;
            plotSingularValues_1Dimensional(arr_SingularValues, limits_k, limits_t);
            
            
            
            
            if strcmp(str_variable , 'x')
    
                % t2 is to be determined
                txtTitle = "Minimum singular values of S_{}";
                txtXLabel = "$k_{2}$";
                txtYLabel = "$\log_{10} \left( \sigma_{t_{1},k_{2}} \right)$";
                
            else
                % t1 is to be determined
                txtTitle = "Minimum singular values of S_{k1,t2}";
                txtXLabel = "$k_{1}$";
                txtYLabel = "$\log_{10}\left(\sigma_{k_{1},t_{2}}\right)$";
            end
            
            
            plotMinimumSingularValues_1Dimensional(vMinimumSingularValues, ...
                limits_k, limits_t, ...
                txtTitle, txtXLabel, txtYLabel);
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
    alpha = 1;
    return;
end


if strcmp(str_variable, 'x')
   
    t2 = t_star;
    
else
    
    t1 = t_star;
    
end




LineBreakLarge();
fprintf('t1 : %i \n', t1 );
fprintf('t2 : %i \n', t2 );
LineBreakLarge();

    
    
    
    
    
    
    
    
end