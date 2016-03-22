function [t, th1, th2] = Get_t(fxy_matrix,gxy_matrix,m,n)
% Get the total degree t of the two input polynomials f(x,y) and g(x,y)

% Initialise the global variables
global BOOL_PREPROC
global BOOL_Q
global PLOT_GRAPHS

% Get degrees of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get degrees of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% Degree elevate f(x,y) so that its coefficients are contained in a square
fxy_matrix_delv = DegreeElevate(fxy_matrix,m-m1,m-m2);

% Degree elevate g(x,y) so that its coefficients are contained in a square
% matrix.
gxy_matrix_delv = DegreeElevate(gxy_matrix,n-n1,n-n2);

% Get degrees of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix_delv);

% Get degrees of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix_delv);

% Initialise some vectors
minmn = min(m,n);

vMinimumSingularVal = zeros(1,minmn);
vCondition = zeros(1,minmn);
vOptTheta1 = zeros(1,minmn);
vOptTheta2 = zeros(1,minmn);
Data_RowNorm = [];
Data_DiagNorm = [];


% for each possible total degree
for k=1:1:min(m,n)
    
    % Apply preprocessing
    switch BOOL_PREPROC
        case 'y'
            
            % Preproecessor One - Normalise by geometric mean
            [lambda, mu] = GetGeometricMean(fxy_matrix_delv,gxy_matrix_delv,k,k);
            
            % Normalise f(x,y)
            fxy_matrix_n = fxy_matrix_delv./lambda;
            
            % Normalise g(x,y)
            gxy_matrix_n = gxy_matrix_delv./mu;
            
            % Preprocessor Two and Three - LinProg to obtain optimal values
            % of alpha, theta_1 and theta_2
            
            % Get the maximum and minimum entries of f(x,y) in the
            % Sylvester matrix S(f,g)
            [max_mtrx_f, min_mtrx_f] = GetMaxMin(fxy_matrix_n,n1,n2,k,k);
            
            % Get the maximum and minimum entries of g(x,y) in the
            % Sylvester matrix S(f,g)
            [max_mtrx_g, min_mtrx_g] = GetMaxMin(gxy_matrix_n,m1,m2,k,k);
            
            % Get optimal values of alpha and theta
            [alpha, th1, th2] = OptimalAlphaTheta(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g);
            
            
        case 'n'
            
            fxy_matrix_n = fxy_matrix_delv;
            gxy_matrix_n = gxy_matrix_delv;
            
            alpha = 1;
            th1 = 1;
            th2 = 1;
            lambda = 1;
            mu = 1;
        otherwise
            error('err')
    end
    
    vOptTheta1(k)       = th1;
    vOptTheta2(k)       = th2;
    opt_alpha_vec(k)    = alpha;
    lambda_vec(k)       = lambda;
    mu_vec(k)           = mu;
    
    
    %% Build the Sylvester matrix S_{k,k}
    
    fww_matrix = GetWithThetas(fxy_matrix_n,th1,th2);
    gww_matrix = GetWithThetas(gxy_matrix_n,th1,th2);
    
    % Build two Cauchy matrices, the first for coefficients of fxy and the
    % second for the coefficients of gxy
    T_f = BuildT1(fww_matrix,n1-k,n2-k);
    T_g = BuildT1(gww_matrix,m1-k,m2-k);
    
    % Build the diagonal matrix D^{-1}
    D = BuildD(k,k,m1,m2,n1,n2);
    
    
    % Include Q / Exclude Q from Sylvester Matrix
    switch BOOL_Q
        case 'y'
            
            % Build the diagonal matrix Q such that Q * [v \\ u] gives the
            % coefficients of u and v in the scaled bernstein basis
            Q1 = BuildQ1(n1-k,n2-k);
            Q2 = BuildQ1(m1-k,m2-k);
                        
            % Build the Two Partitions of the Syvlester Matrix 
            % S = [C(f) C(g)].
            C_f = D * T_f * Q1;
            C_g = D * T_g * Q2;
        case 'n'
            % Do Nothing
        otherwise
            error('BOOL_Q is either (y) or (n)')
    end
    
    % Build the Sylvester matrix
    Sk = [C_f  alpha.* C_g];
    
    ratio_max_min_entries_proc(k) = max(max(Sk))./min(min(Sk));
    
    
    %% Get QR Decomposition
    % Using QR Decomposition of the sylvester matrix
    [~,R] = qr(Sk);
    
    % Take absolute values.
    R = abs(R);
    
    % Get number of rows in R1
    [R1_rows,~] = size(diag(R));
    
    % Obtain R1 the top square of the R matrix.
    R1 = R(1:R1_rows,1:R1_rows);
    
    % Get Norms of each row in the matrix R1
    R1_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get ONLY the diagonal elements and normalise them.
    R1_DiagNorm = diag(R1)./norm(diag(R1));
    
    % Build Array of all diagonals
    % Scatter Plot Data
    ks = k.*ones(size(R1_RowNorm));
    ns = 1:1:size(R1_RowNorm,1);
    
    X = [ks R1_RowNorm ns'];
    Data_RowNorm = [Data_RowNorm; X];
    
    X2 = [ks R1_DiagNorm ns'];
    Data_DiagNorm = [Data_DiagNorm;X2];
    
    %%
    
    % Get SVD of unproc and processed Sylvester Surbesultant S_{k,k}
    vMinimumSingularVal(k) = min(svd(Sk));
    %min_sing_val_vec_unproc(k) = min(svd(Sk_unproc));
    
    % Get the condition of Sk
    vCondition(k) = cond(Sk);
    %cond_vec_unproc(k) = cond(Sk_unproc);
    
end


switch PLOT_GRAPHS
    case 'y'
        %% plot the minimum singular values
        figure('name','Min Sing Val')
        title('minimum Singular Value for each subresultant matrix S_{k,k}')
        hold on
        plot(log10(vMinimumSingularVal),'-s','DisplayName','Preprocessed');
        %plot(log10(min_sing_val_vec_unproc),'-s','DisplayName','Unprocessed');
        xlabel('k : index of subresultant')
        legend(gca,'show')
        ylabel('log_{10} Minimum Singular Value')
        
        hold off
        %% Plot the condition numbers of each subresultant
        figure('name','Condition S_{k}')
        title('Condition of each subresultant S_{k,k}')
        hold on
        plot(log10(vCondition),'-s','DisplayName','Preprocessed');
        %plot(log10(cond_vec_unproc),'-s','DisplayName','Unprocessed');
        xlabel('k : index of subresultant S_{k}')
        ylabel('log_{10} Condition Number')
        legend(gca,'show')
        hold off
        
        %%
        figure('name','Diag Norms')
        plot(Data_DiagNorm(:,1),(log10(Data_DiagNorm(:,2))),'*')
        axis([0.9,min(m,n),-inf,+inf])
        xlabel('k')
        ylabel('Normalised Row Sums of R1 in S_{k}')
        title(['Normalised Row Sums of R1 fom the QR decomposition of each subresultant S_{k} \newline '...
            'm = ' int2str(m) ', n = ' int2str(n) '(Original)']);
        hold off
    case 'n'
    otherwise
        error('error: plot_graphs must be either (y) or (n)')
end
%%
[svd_val,svd_maxindex] = max(diff(log10(vMinimumSingularVal)));

% [rowdiag_val,rowdiag_maxindex] = min(diff(log10(ratio_maxmin_diag_vec)));
% fprintf('Total Degree Calculated By Max:Min Row Diags: %i \n',rowdiag_maxindex);
%
% [rowsum_val,rowsum_maxindex] = min(diff(log10(ratio_maxmin_rowsum_vec)));
% fprintf('Total Degree Calculated By Max:Min Row Sums: %i \n',rowsum_maxindex);

fprintf('-----------------------------------------------------------------\n')


val = svd_val;
index = svd_maxindex;

% check if the maximum change is significant
threshold = 3.5;
if abs(val) < threshold
    % not significant
    fprintf('Differences between minimum singular values of S_{k,k} are not significant \n')
    t = min(m,n);
else
    % change is significant
    t = index;
    fprintf('Total Degree Calculated By Minimum Singular Values: %i \n',t);
    
end


% Set the optimal theta 1 and theta 2
th1 = vOptTheta1(t);
th2 = vOptTheta2(t);

fprintf('---------------------------------------------------------\n')
fprintf('\n')
fprintf('Total Degree of GCD as calculated by GetDegree_Total() : \n')
fprintf('Tota Degree = %i',t)
fprintf('\n')
fprintf('---------------------------------------------------------\n')


end