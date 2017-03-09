function Q1 = BuildQ1_Univariate(n_k)
% BuildQ1_Univariate(n_k)
%
% Build the partition of the matrix Q, corresponding to the binomial
% coefficients of v(y). (Note this function works for both v(y) and u(y)
%
% % Inputs
%
% n_k : (Int) where n_k = n-k
%
% n : (Int) Degree of polynomial g(y)
%
% k : (Int) Index of subresultant S_{t} to be built.

% Produce the matrix Q1
Q1 = diag(GetBinomials(n_k));

end