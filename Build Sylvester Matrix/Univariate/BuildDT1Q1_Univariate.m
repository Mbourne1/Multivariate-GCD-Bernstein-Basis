function DT1Q1 = BuildDT1Q1_Univariate(fx, n_k)
% BuildDT1Q1_Univariate(fx, n_k)
%
% Build Toeplitz matrix for Sylvester Matrix Partitions
%
%
% Inputs.
%
% fx  :  Coefficients of polynomial f given in Bernstein basis n  : degree
%          of polynomial g k  : the index of the subresultant being built
%
% n_k :  Degree of polynomial which f is to be multiplied by.


% Get the degree of polynomial f(x)
m = GetDegree_Univariate(fx);

% Build matrices D^{-1}
D = BuildD_Univariate_2Polys(m, n_k);

% Build the matrix T1
T1 = BuildT1_Univariate(fx, n_k);

% Build the matrix Q1
Q1 = BuildQ1_Univariate(n_k);

% Get the matrix D^{-1} * T_{n-k}(f) * Q_{n-k}
DT1Q1 = D*T1*Q1;
        
        
end










