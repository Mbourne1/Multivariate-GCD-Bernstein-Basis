function [t1_possible, t2_possible, lambda, mu, rho, th1, th2] = GetGCDDegree_Bivariate_3Polys_WithDegreeElevation(fxy, gxy, hxy)
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


% Get degree of f(x,y) and g(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);




% Set my_limits : Always compute all Sylvester subresultant matrices. Use
% limits_t as an indicator of where the GCD degree should be.
method = 1;

switch method 
    case 1
         
        m_star = max(m1, m2);
        n_star = max(n1, n2);
        o_star = max(o1, o2);
        
        %m_star = m1 + m2;
        %n_star = n1 + n2;
        
        p1 = m_star - m1;
        p2 = m_star - m2;
        q1 = n_star - n1; 
        q2 = n_star - n2;
        r1 = o_star - o1;
        r2 = o_star - o2;
        
        fprintf('Degree Elevations f wrt x %i \n', p1);
        fprintf('Degree Elevations f wrt y %i \n', p2);
        fprintf('Degree Elevations g wrt x %i \n', q1);
        fprintf('Degree Elevations g wrt y %i \n', q2);
        fprintf('Degree Elevations h wrt x %i \n', r1);
        fprintf('Degree Elevations h wrt y %i \n', r2);
        
        limits_k = [1 min([m_star, n_star, o_star])];
        
    case 2
        p1 = 0;
        p2 = 0;
        q1 = 0; 
        q2 = 0;
        r1 = 0;
        r2 = 0;
        k_star = min(m1,n1) + min(m2,n2) ;
        
        limits_k = [0 k_star];
        
end


% Degree elevate f(x,y) so that its coefficients are contained in a square
fxy_matrix_delv = DegreeElevate_Bivariate(fxy, p1, p2);
gxy_matrix_delv = DegreeElevate_Bivariate(gxy, q1, q2);
hxy_matrix_delv = DegreeElevate_Bivariate(hxy, r1, r2);
        
limits_t = limits_k;

lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% Set the number of subresultants to be built
nSubresultants = upperLimit_k - lowerLimit_k + 1;

% if upper limit is equal to lower limit
if upperLimit_k == lowerLimit_k
    t_star = upperLimit_k;
    th1 = 1;
    th2 = 1;
    return;
end



% Initialise some vectors.
vCondition = zeros(nSubresultants, 1);

% Initialise vectors
vGM_fx = zeros(nSubresultants, 1);
vGM_gx = zeros(nSubresultants, 1);

vLambda = zeros(nSubresultants, 1);
vMu = zeros(nSubresultants, 1);
vRho = zeros(nSubresultants, 1);

vTh1 = zeros(nSubresultants, 1);
vTh2 = zeros(nSubresultants, 1);

% Initialise cell arrays
arr_Sk = cell(nSubresultants, 1);
arr_R1 = cell(nSubresultants, 1);
arr_SingularValues = cell(nSubresultants, 1);

% for each possible total degree
for i = 1 : 1 : nSubresultants
    
    k = lowerLimit_k + (i-1);
    
    % Apply preprocessing
    [GM_fx, GM_gx, GM_hx, lambda, mu, rho, th1, th2] = ...
        Preprocess_Bivariate_3Polys(fxy_matrix_delv, gxy_matrix_delv, hxy_matrix_delv, k, k);
    
    vGM_fx(i) = GM_fx;
    vGM_gx(i) = GM_gx;
    vGM_hx(i) = GM_hx;
    vLambda(i) = lambda;
    vMu(i) = mu;
    vRho(i)= rho;
    
    vTh1(i) = th1;
    vTh2(i) = th2;
    
    % Get f(x,y) normalized by geometric mean lambda.
    fxy_matrix_n = fxy_matrix_delv ./ vGM_fx(i);
    gxy_matrix_n = gxy_matrix_delv ./ vGM_gx(i);
    hxy_matrix_n = hxy_matrix_delv ./ vGM_hx(i);
    
    % Get f(w,w) from f(x,y)
    % Get g(w,w) from g(x,y)
    fww_matrix = GetWithThetas(fxy_matrix_n, vTh1(i), vTh2(i));
    gww_matrix = GetWithThetas(gxy_matrix_n, vTh1(i), vTh2(i));
    hww_matrix = GetWithThetas(hxy_matrix_n, vTh1(i), vTh2(i));
    
    % Build the k-th subresultant matrix.
    arr_Sk{i} = BuildSubresultant_Bivariate_3Polys(lambda.*fww_matrix, mu.*gww_matrix, rho.*hww_matrix, k, k);
    
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
            
            txtXlabel = "k";
            txtYlabel = "\sigma_{k,k}";
            txtTitle = "Singular Values of Sk";
            
            plotMinimumSingularValues_1Dimensional(vMinimumSingularValues,...
                limits_k, limits_t, txtTitle, txtXlabel, txtYlabel);
            
            
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
    mu = 1;
    rho = 1;
    return;
end



% Set the optimal theta 1 and theta 2
th1 = vTh1(t_star - lowerLimit_k + 1);
th2 = vTh2(t_star - lowerLimit_k + 1);
lambda = vLambda(t_star - lowerLimit_k + 1);
mu = vMu(t_star - lowerLimit_k + 1);
rho = vRho(t_star - lowerLimit_k + 1);

fprintf('%s : Value of t* is given by : %i \n', mfilename, t_star)

t1_possible = t_star - (min([p1, q1, r1]));
t2_possible = t_star - (min([p2, q2, r2]));

LineBreakLarge();
fprintf('Either \n')
fprintf('t1 : %i \n', t1_possible );
fprintf('t2 : %i \n', t2_possible );
LineBreakLarge();





end

