function [t, th1, th2] = GetGCDDegree_Total(fxy_matrix,gxy_matrix,m,n)
% Get the total degree t of the two input polynomials f(x,y) and g(x,y)


lower_lim = 1;
upper_lim = min(m,n);

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
vMinimumSingularValues = zeros(min_mn,1);
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
        Preprocess(fxy_matrix_delv,gxy_matrix_delv,k,k);
    
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
    
    vSingularValues = svd(Sk);
    
    % Get SVD of unproc and processed Sylvester Surbesultant S_{k,k}
    vMinimumSingularValues(k) = min(vSingularValues);
    
    % Get the condition of Sk
    vCondition(k) = cond(Sk);
    
end

% Check if more than one subresultant exists
if min(m,n) == 1
    t = GetGCDDegree_OneSubresultant(vSingularValues);
else
    t = GetGCDDegree_MultipleSubresultants(vMinimumSingularValues,[lower_lim, upper_lim]);
    PlotGraphs_TotalDegree() 
end

if t == 0
    th1 = 1;
    th2 = 1;
    return;
end
%%


% Set the optimal theta 1 and theta 2
th1 = vTh1(lower_lim + t - 1);
th2 = vTh2(lower_lim + t - 1);

LineBreakLarge();
fprintf([mfilename ' : ' sprintf('Total degree : %i','t')])
LineBreakLarge()

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

