function [t1,t2,lambda,mu,alpha,th1,th2] = GetGCDDegree_Relative(fxy,gxy)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y)
%
% % Inputs.
%
%   fxy : Coefficient matrix of polynomial f(x,y)
%
%   gxy : Coefficient matrix of polynomial g(x,y)
%
% % Outputs 
%
%   t1 : Degree of d(x,y) with respect to x
% 
%   t2 : Degree of d(x,y) with respect to y
%
%   lambda : Geometric mean of entries in first partition of Syvlester
%   matrix
%
%   mu : Geometric mean of Entries in second partition of Sylvester matrix
%
%   alpha : Optimal value of alpha
%
%   th1 : Optimal value of theta_{1}
%
%   th2 : Optimal value of theta_{2}


%
global SETTINGS


% Get the degree structure of polynomial f(x,y)
[m1,m2] = GetDegree(fxy);

% Get the degree structure of polynomial g(x,y)
[n1,n2] = GetDegree(gxy);

% Get all k1 k2 pairs
k1k2Pairs = GetPairs_All(m1,m2,n1,n2);

% If only one entry in matrix of possible (k1,k2) pairs exist
% then (t1,t2) = (k1,k2)
[nPairs,~] = size(k1k2Pairs);

% if nPairs == 1
%     fprintf('Only one possible combination of (t1,t2) \n')
%     t1 = k1k2Pairs(1,1);
%     t2 = k1k2Pairs(1,2);
%     return
% end

% Initialise some vectors
vMinimumSingularValues_all = zeros(nPairs,1);
vAlpha_all = zeros(nPairs,1);
vTheta1_all = zeros(nPairs,1);
vTheta2_all = zeros(nPairs,1);
vLambda_all = zeros(nPairs,1);
vMu_all = zeros(nPairs,1);

% For every row in the matrix
for i = 1:1:nPairs
    
    % Get the ith pair of k_{1} and k_{2} values
    v_k1 = k1k2Pairs(i,1);
    v_k2 = k1k2Pairs(i,2);
    
    % Preprocessing
    [lambda, mu, alpha, th1,th2] = Preprocess(fxy,gxy,v_k1,v_k2);
    
    % Divide f(x) by geometric mean
    fxy_matrix_n = fxy ./lambda;
    
    % Divide g(x) by geometric mean
    gxy_matrix_n = gxy ./mu;
    
    % Get f(w,w) from f(x,y)
    fww_matrix = GetWithThetas(fxy_matrix_n,th1,th2);
    
    % Get g(w,w) from g(x,y)
    gww_matrix = GetWithThetas(gxy_matrix_n,th1,th2);
    
    % Build the k1,k2 subresultant matrix
    Sk1k2 = BuildDTQ(fww_matrix,alpha.*gww_matrix,v_k1,v_k2);
    
    % Get the singular values
    vSingularValues = svd(Sk1k2);
    
    % Store the minimum Singular value
    vMinimumSingularValues_all(i) = min(vSingularValues);
    % Store the optimal alpha
    vAlpha_all(i) = alpha;
    % Store the optimal theta_{1}
    vTheta1_all(i) = th1;
    % Store the optimal theta_{2}
    vTheta2_all(i) = th2;
    % Store the optimal lambda
    vLambda_all(i) = lambda;
    % Store the optimal mu
    vMu_all(i) = mu;
    
end


% Sort all k1 k2 pairs by their total k_{total} =  k1 + k2. Get the minimum of the
% minimum singular values for each total.

% Get the sum of k1 + k2 for all (k1,k2) pairs
sumk1k2 = sum(k1k2Pairs,2);

% Create a matrix of data
data = [k1k2Pairs, sumk1k2, vMinimumSingularValues_all,...
    vAlpha_all, vTheta1_all, vTheta2_all, vLambda_all, vMu_all];

% Get the minimum sum
min_val = min(sumk1k2);

% Get the maximum sum
max_val = max(sumk1k2);

% Get the number of possible sum values.
nValues = max_val - min_val + 1;

% Initialise a vector of k1
v_k1 = zeros(nValues,1);

% Initialise a vector of k2
v_k2 = zeros(nValues,1);

% Initialise a vector to store minimum singular values
vMinimumSingularValues = zeros(nValues,1);


vAlpha = zeros(nValues,1);
vTheta1 = zeros(nValues,1);
vTheta2 = zeros(nValues,1);
vLambda = zeros(nValues,1);
vMu = zeros(nValues,1);



for i = min_val:1:max_val
    
    k = i - min_val + 1;
    
    data_filtered = data(data(:, 3) == i, :);
    
    % Get the minimum singular value of all the [k1,k2] pairs
    [~,index] = min(data_filtered(:,4));
    
    v_k1(k) = data_filtered(index,1);
    v_k2(k) = data_filtered(index,2);
    vMinimumSingularValues(k) = data_filtered(index,4);
    vAlpha(k) = data_filtered(index,5);
    vTheta1(k) = data_filtered(index,6);
    vTheta2(k) = data_filtered(index,7);
    vLambda(k) = data_filtered(index,8);
    vMu(k) = data_filtered(index,9);
    
end

% Plot the minimum singular values
PlotGraphs_DegreeRelative();


% If only one value (t1+t2)
if (nValues ==1)
    t1 = v_k1;
    t2 = v_k2;
    alpha = vAlpha(1);
    th1 = vTheta1(1);
    th2 = vTheta2(1);
    lambda = vLambda(1);
    mu = vMu(1);
    return
end



% Get maximum change in Singular values
[maxChange,index] = max(diff(log10(vMinimumSingularValues)));

if (maxChange < SETTINGS.THRESHOLD_RANK)
    fprintf([mfilename ' : ' 'Insignificant Change \n'])
    fprintf([mfilename ' : ' 'All subresultants are either full rank or rank deficient \n'])
    fprintf([mfilename ' : ' 'All subresultants are full rank \n'])
    t1 = v_k1(end);
    t2 = v_k2(end);
    alpha = vAlpha(end);
    th1 = vTheta1(end);
    th2 = vTheta2(end);
    lambda = vLambda(end);
    mu = vMu(end);
    
else
    fprintf([mfilename ' : ' 'Significant Change \n'])
    t1 = v_k1(index);
    t2 = v_k2(index);
    alpha = vAlpha(index);
    th1 = vTheta1(index);
    th2 = vTheta2(index);
    lambda = vLambda(index);
    mu = vMu(index);
end


LineBreakMedium()
fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n',t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n',t2)])
LineBreakMedium()

end




