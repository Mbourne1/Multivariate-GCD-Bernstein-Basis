function [t, th1, th2] = GetGCDDegree_Total2(fxy_matrix,gxy_matrix,m,n, limits_t)
% Get the total degree t of the two input polynomials f(x,y) and g(x,y)

% Get degrees of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get degrees of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% Set the lower limit for the value of t
lower_limit = limits_t(1);

% Set the upper limit for the value of t
upper_limit = limits_t(2);

% Set the number of subresultants to be built
n_subresultants = upper_limit - lower_limit +1;

% if upper limit is equal to lower limit
if upper_limit == lower_limit
   t = upper_limit;
   th1 = 1;
   th2 = 1;
   return; 
end

% Degree elevate f(x,y) so that its coefficients are contained in a square
fxy_matrix_delv = DegreeElevate(fxy_matrix,m-m1,m-m2);

% Degree elevate g(x,y) so that its coefficients are contained in a square
% matrix.
gxy_matrix_delv = DegreeElevate(gxy_matrix,n-n1,n-n2);

% Initialise some vectors.
vMinimumSingularValues = zeros(n_subresultants,1);
vCondition = zeros(n_subresultants,1);

vGM_fx = zeros(n_subresultants,1);
vGM_gx = zeros(n_subresultants,1);
vAlpha = zeros(n_subresultants,1);
vTh1 = zeros(n_subresultants,1);
vTh2 = zeros(n_subresultants,1);


Data_RowNorm = [];
Data_DiagNorm = [];




% for each possible total degree
for k = lower_limit : 1 : upper_limit
    
    i = k - lower_limit + 1;
    
    % Apply preprocessing
    [vGM_fx(i),vGM_gx(i),vAlpha(i), vTh1(i),vTh2(i)] = ...
        Preprocess(fxy_matrix_delv,gxy_matrix_delv,k,k);
    
    % Get f(x,y) normalized by geometric mean lambda.
    fxy_matrix_n = fxy_matrix_delv./vGM_fx(i);
    
    % Get g(x,y) normalised by geometric mean mu.
    gxy_matrix_n = gxy_matrix_delv./vGM_gx(i);
    
    % Get f(w,w) from f(x,y)
    fww_matrix = GetWithThetas(fxy_matrix_n,vTh1(i),vTh2(i));
    
    % Get g(w,w) from g(x,y)
    gww_matrix = GetWithThetas(gxy_matrix_n,vTh1(i),vTh2(i));
    
    % Build the k-th subresultant matrix.
    Sk = BuildDTQ(fww_matrix,vAlpha(i).*gww_matrix,k,k);
    
    % Get the matrix R1 from QR decomposition of S_{k}
    R1 = GetR1(Sk);
    
    % Get Norms of each row in the matrix R1
    R1_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get ONLY the diagonal elements and normalise them.
    R1_DiagNorm = diag(R1)./norm(diag(R1));
       
    Data_RowNorm = AddToData(R1_RowNorm,Data_RowNorm,k);
    Data_DiagNorm = AddToData(R1_DiagNorm,Data_DiagNorm,k);
    
    % Get SVD of unproc and processed Sylvester Surbesultant S_{k,k}
    vMinimumSingularValues(k) = min(svd(Sk));
    
    % Get the condition of Sk
    vCondition(k) = cond(Sk);
    
end

PlotGraphs()


[svd_val,svd_maxindex] = max(diff(log10(vMinimumSingularValues)));


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
