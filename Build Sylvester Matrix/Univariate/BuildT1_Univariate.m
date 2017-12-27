function T1 = BuildT1_Univariate(fx,n_k)
% BuildT1_Univariate(fx,n_k)
%
% Build a Toeplitz Matrix of coefficients of f(x).
% T1 \in \mathbb{R}^{(m+n-k+1)\times(n-k+1)}
%
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% n_k : (Int) Degree of polynomial v_{k}(x)
%
% n : (Int) Degree of polynomial g(x)
%
% k : (Int) Index of subresultant S_{k} to be formed. 
%
% % Outputs
%
% T1 : (Matrix) Convolution matrix T_{n-k}(f(x))



% Get degree of polynomail f(x)
m = GetDegree_Univariate(fx);

% Get f(w) with binomial coefficients;
fx_bi = GetWithBinomials_Univariate(fx);

% Initialise empty matrix T1, for storing Toeplitz T_{k}(f)
T1 = zeros(m+n_k+1, n_k+1);

% for each column of T1
for j = 0:1:n_k
    
    T1(j + 1 : m + j + 1, j + 1) = fx_bi;
    
end


end