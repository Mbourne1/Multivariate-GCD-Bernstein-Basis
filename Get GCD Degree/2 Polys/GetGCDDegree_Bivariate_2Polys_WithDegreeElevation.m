function [t_star, th1, th2] = GetGCDDegree_Bivariate_2Polys_WithDegreeElevation(fxy, gxy, p1, p2, q1, q2, limits_t)
% GetGCDDegree_Total_Bivariate_2Polys(fxy, gxy, m, n, limits_t)
%
% Get the total degree t of the two input polynomials f(x,y) and g(x,y)
%
% % Inputs
%
% fxy : Coefficients of the polynomial f(x,y)
%
% gxy : Coefficients of the polynomial g(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
%
% limits_t :
%
% % Outputs
%
% t :
%
% th1 :
%
% th2 :

global SETTINGS



% Degree elevate f(x,y) so that its coefficients are contained in a square
fxy_matrix_delv = DegreeElevate_Bivariate(fxy, p1, p2);


% Degree elevate g(x,y) so that its coefficients are contained in a square
% matrix.
gxy_matrix_delv = DegreeElevate_Bivariate(gxy, q1, q2);



% Set my_limits : Always compute all Sylvester subresultant matrices. Use
% limits_t as an indicator of where the GCD degree should be.
myLimits = [0 min(m_star, n_star)];

myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

% Set the number of subresultants to be built
nSubresultants = myUpperLimit - myLowerLimit + 1;

% if upper limit is equal to lower limit
if myUpperLimit == myLowerLimit
    t_star = myUpperLimit;
    th1 = 1;
    th2 = 1;
    return;
end



% Initialise some vectors.
vCondition = zeros(nSubresultants, 1);

% Initialise vectors
vGM_fx = zeros(nSubresultants, 1);
vGM_gx = zeros(nSubresultants, 1);
vAlpha = zeros(nSubresultants, 1);
vTh1 = zeros(nSubresultants, 1);
vTh2 = zeros(nSubresultants, 1);

% Initialise cell arrays
arr_Sk = cell(nSubresultants, 1);
arr_R1 = cell(nSubresultants, 1);
arr_SingularValues = cell(nSubresultants, 1);

% for each possible total degree
for i = 1:1:nSubresultants
    
    k = myLowerLimit + (i-1);
    
    % Apply preprocessing
    [vGM_fx(i),vGM_gx(i),vAlpha(i), vTh1(i),vTh2(i)] = ...
        Preprocess(fxy_matrix_delv, gxy_matrix_delv, k, k);
    
    % Get f(x,y) normalized by geometric mean lambda.
    fxy_matrix_n = fxy_matrix_delv./vGM_fx(i);
    
    % Get g(x,y) normalised by geometric mean mu.
    gxy_matrix_n = gxy_matrix_delv./vGM_gx(i);
    
    % Get f(w,w) from f(x,y)
    fww_matrix = GetWithThetas(fxy_matrix_n, vTh1(i), vTh2(i));
    
    % Get g(w,w) from g(x,y)
    gww_matrix = GetWithThetas(gxy_matrix_n, vTh1(i), vTh2(i));
    
    % Build the k-th subresultant matrix.
    arr_Sk{i} = BuildDTQ_Bivariate_2Polys_NewMethod(fww_matrix, vAlpha(i).*gww_matrix, 1, k);
    
    % Get the matrix R1 from QR decomposition of S_{k}
    arr_R1{i} = GetR1(arr_Sk{i});
    
    % Get array of singular values
    arr_SingularValues{i} = svd(arr_Sk{i});
    
    % Get the condition of Sk
    vCondition(i) = cond(arr_Sk{i});
    
end

fprintf('Metric used to compute the degree of the GCD : %s \n', SETTINGS.RANK_REVEALING_METRIC)
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals
switch SETTINGS.RANK_REVEALING_METRIC
    
    case 'Singular Values'
        
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
            %plotSingularValues_TotalDegree(arr_SingularValues);
            plotMinimumSingularValues_TotalDegree(vMinimumSingularValues, myLimits);
        end
        
    case 'R1 Row Norms'
        
        % Initialise vectors to store Max and Min Row Norms of R1_{k}
        arr_R1_RowNorm = cell(nSubresultants,1);
        vMaxRowNorm = zeros(nSubresultants,1);
        vMinRowNorm = zeros(nSubresultants,1);
        
        
        for i = 1:1:nSubresultants
            
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
            plotMaxMinDiagonals_R1_Total(vMaxDiagonal_R1, vMinDiagonal_R1, myLimits, limits_t);
            %plotDiagonals_R1(arr_R1_diag);
        end
        
    case 'Residuals'
        error('Code not complete')
    otherwise
        error('err');
        
end



if myUpperLimit == myLowerLimit
    
    t_star = GetGCDDegree_OneSubresultant(metric, myLimits);
    
    
    
else
    
    t_star = GetGCDDegree_MultipleSubresultants(metric, myLimits);
    
    
    
end

if t_star == 0
    th1 = 1;
    th2 = 1;
    return;
end
%%


% Set the optimal theta 1 and theta 2
th1 = vTh1(t_star - myLowerLimit + 1);
th2 = vTh2(t_star - myLowerLimit + 1);

LineBreakLarge();
fprintf('%s : Value of t* is given by : %i \n', mfilename, t_star)
LineBreakLarge();

end

