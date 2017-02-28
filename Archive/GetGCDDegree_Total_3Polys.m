function [t, th1, th2] = GetGCDDegree_Total_3Polys(fxy, gxy, hxy, m, n, o, limits_t)
% Get the total degree t of the two input polynomials f(x,y) and g(x,y)
%
% % Inputs
%
% [fxy, gxy, hxy] : Coefficients of the polynomial f(x,y) and g(x,y)
%
% [m, n, o] : Total degree of f(x,y) and g(x,y)
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



% Get degrees of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degrees of polynomial g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Get degrees of polynomial h(x,y)
[o1, o2] = GetDegree_Bivariate(hxy);

bool_method = '1';
switch bool_method 
    case '1'
        m = max(m1, m2);
        n = max(n1, n2);
        o = max(o1, o2);
        upper_lim = min(m,n,o);
    case '2'
        %m = m;
        %n = n;
        upper_lim = limits_t(2);
end

% Set the lower limit for the value of t
lower_lim = limits_t(1);


% Set the number of subresultants to be built
n_subresultants = upper_lim - lower_lim +1;

% if upper limit is equal to lower limit
if upper_lim == lower_lim
   t = upper_lim;
   th1 = 1;
   th2 = 1;
   return; 
end

% Degree elevate f(x,y) so that its coefficients are contained in a square
fxy_matrix_delv = DegreeElevate(fxy, m-m1, m-m2);

% Degree elevate g(x,y) so that its coefficients are contained in a square
% matrix.
gxy_matrix_delv = DegreeElevate(gxy, n-n1, n-n2);


hxy_matrix_delv = DegreeElevate(hxy, o-o1, o-o2);

% Initialise some vectors.
vMinimumSingularValues = zeros(n_subresultants,1);
vCondition = zeros(n_subresultants,1);

vGM_fx = zeros(n_subresultants,1);
vGM_gx = zeros(n_subresultants,1);
vGM_hx = zeros(n_subresultants,1);

vAlpha = zeros(n_subresultants,1);
vTh1 = zeros(n_subresultants,1);
vTh2 = zeros(n_subresultants,1);


Data_RowNorm = [];
Data_DiagNorm = [];




% for each possible total degree
for k = lower_lim : 1 : upper_lim
    
    i = k - lower_lim + 1;
    
    % Apply preprocessing
    %[vGM_fx(i),vGM_gx(i),vAlpha(i), vTh1(i),vTh2(i)] = ...
    %    Preprocess(fxy_matrix_delv,gxy_matrix_delv,k,k);
    
    % Temporary Measure
    fprintf('')
    vGM_fx(i) = 1;
    vGM_gx(i) = 1;
    vGM_hx(i) = 1;
    vAlpha(i) = 1;
    vTh1(i) = 1;
    vTh2(i) = 1;
    
    % Get f(x,y) normalized by geometric mean lambda.
    fxy_matrix_n = fxy_matrix_delv./vGM_fx(i);
    
    % Get g(x,y) normalised by geometric mean mu.
    gxy_matrix_n = gxy_matrix_delv./vGM_gx(i);
    
    % Get h(x,y) normalised by geometric mean
    hxy_matrix_n = hxy_matrix_delv./vGM_hx(i);
    
    
    % Get f(w,w) from f(x,y)
    fww_matrix = GetWithThetas(fxy_matrix_n, vTh1(i), vTh2(i));
    
    % Get g(w,w) from g(x,y)
    gww_matrix = GetWithThetas(gxy_matrix_n, vTh1(i), vTh2(i));
    
    % Get h(\omega_{1},\omega_{2}) from h(x,y)
    hww_matrix = GetWithThetas(hxy_matrix_n, vTh1(i), vTh2(i));
    
    % Build the k-th subresultant matrix.
    Sk = BuildDTQ_3Polys(fww_matrix,vAlpha(i).*gww_matrix, hww_matrix, k, k);
    
    % Get the matrix R1 from QR decomposition of S_{k}
    R1 = GetR1(Sk);
    
    % Get Norms of each row in the matrix R1
    R1_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get ONLY the diagonal elements and normalise them.
    R1_DiagNorm = diag(R1)./norm(diag(R1));
       
    Data_RowNorm = AddToData(R1_RowNorm,Data_RowNorm,k);
    Data_DiagNorm = AddToData(R1_DiagNorm,Data_DiagNorm,k);
    
    % Get SVD of unproc and processed Sylvester Surbesultant S_{k,k}
    vMinimumSingularValues(i) = min(svd(Sk));
    
    % Get the condition of Sk
    vCondition(i) = cond(Sk);
    
end

if min(m,n) == 1
    t = GetGCDDegree_OneSubresultant(vSingularValues);
    
    
else
    PlotGraphs_TotalDegree2()
    t = GetGCDDegree_MultipleSubresultants(vMinimumSingularValues,[lower_lim,upper_lim]);
    
    
end

if t == 0
    th1 = 1;
    th2 = 1;
    return;
end
%%


% Set the optimal theta 1 and theta 2
th1 = vTh1(t - lower_lim + 1);
th2 = vTh2(t - lower_lim + 1);


fprintf('%s : Total Degree of GCD : %i \n',mfilename,t)


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
