function [t, th1, th2] = Get_t(fxy_matrix,gxy_matrix,m,n)
% Get the total degree t of the two input polynomials f(x,y) and g(x,y)

% Initialise the global variables
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

% Get min(m,n)
min_mn = min(m,n);

% Initialise some vectors.
vMinimumSingularVal = zeros(min_mn,1);
vCondition = zeros(min_mn,1);

vGM_fx = zeros(min_mn,1);
vGM_gx = zeros(min_mn,1);
vAlpha = zeros(min_mn,1);
vTh1 = zeros(min_mn,1);
vTh2 = zeros(min_mn,1);


Data_RowNorm = [];
Data_DiagNorm = [];


% for each possible total degree
for k=1:1:min_mn
    
    % Apply preprocessing
    
    [vGM_fx(k),vGM_gx(k),vAlpha(k), vTh1(k),vTh2(k)] = ...
        Preprocess(fxy_matrix_delv,gxy_matrix_delv,k);
    
    % Get f(x,y) normalized by geometric mean lambda.
    fxy_matrix_n = fxy_matrix_delv./vGM_fx(k);
    
    % Get g(x,y) normalised by geometric mean mu.
    gxy_matrix_n = gxy_matrix_delv./vGM_gx(k);
    
    % Get f(w,w) from f(x,y)
    fww_matrix = GetWithThetas(fxy_matrix_n,vTh1(k),vTh2(k));
    
    % Get g(w,w) from g(x,y)
    gww_matrix = GetWithThetas(gxy_matrix_n,vTh1(k),vTh2(k));
    
    % Build the k-th subresultant matrix.
    Sk = BuildDTQ(fww_matrix,vAlpha(k).*gww_matrix,k,k);
    
    % Get the matrix R1 from QR decomposition of S_{k}
    R1 = GetR1(Sk);
    
    % Get Norms of each row in the matrix R1
    R1_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get ONLY the diagonal elements and normalise them.
    R1_DiagNorm = diag(R1)./norm(diag(R1));
       
    Data_RowNorm = AddToData(R1_RowNorm,Data_RowNorm,k);
    Data_DiagNorm = AddToData(R1_DiagNorm,Data_DiagNorm,k);
    
    % Get SVD of unproc and processed Sylvester Surbesultant S_{k,k}
    vMinimumSingularVal(k) = min(svd(Sk));
    
    % Get the condition of Sk
    vCondition(k) = cond(Sk);
    
    
    
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
th1 = vTh1(t);
th2 = vTh2(t);

fprintf('---------------------------------------------------------\n')
fprintf('\n')
fprintf('Total Degree of GCD as calculated by GetDegree_Total() : \n')
fprintf('Tota Degree = %i',t)
fprintf('\n')
fprintf('---------------------------------------------------------\n')


end

function [data] = AddToData(vResults,data,k)

% Get number of entries in results vector
[nEntries,~] = size(vResults);

% Get vector of k
vec_k = k.*ones(nEntries,1);

% Get index vector
vec_idx = (1:1:nEntries)';

data = ...
    [
    data;
    vec_k vResults vec_idx;
    ];



end


function R1 = GetR1(Sk)
% Given the k-th Sylvester matrix get the upper triangular square matrix
% R1.
%
% Inputs.
%
% Sk

% Get QR Decomposition
% Using QR Decomposition of the sylvester matrix
[~,R] = qr(Sk);

% Take absolute values.
R = abs(R);

% Get number of rows in R1
[R1_rows,~] = size(diag(R));

% Obtain R1 the top square of the R matrix.
R1 = R(1:R1_rows,1:R1_rows);

end