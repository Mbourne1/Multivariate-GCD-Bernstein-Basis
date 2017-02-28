function T1 = BuildT1_Univariate(fx,n_k)
% BuildT1_univar(fx,n_k)
%
% Build a Toeplitz Matrix of coefficients of f(x).
% T1 \in \mathbb{R}^{(m+n-k+1)\times(n-k+1)}
%
%
% % Inputs
%
% fx : Coefficients of polynomial f(x)
%
% n : Degree of polynomial g
%
% k : Index of subresultant S_{k} to be formed. 


% Get degree of polynomail f(\omega,\theta)
m = GetDegree_Univariate(fx);

% Initialise empty matrix T1, for storing Toeplitz T_{k}(f)
T1 = zeros(m+n_k+1,n_k+1);

% Get f(w) with binomial coefficients;
fx_bi = GetWithBinomials_Univariate(fx);

% for each column of T1
for j = 0:1:n_k
    T1(j+1:m+j+1,j+1) = fx_bi;
end


end