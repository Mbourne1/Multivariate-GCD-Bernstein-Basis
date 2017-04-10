function D = BuildD_Univariate_2Polys(m, n_k)
% Build the matrix D^{-1} which is the diagonal matrix of binomial 
% coefficients. Used to construct the Sylvester Subresultant matrices in 
% format D_{k}^{-1}*T_{k}(f,g)*Q.
% 
%
% % Input
%
% m : (Int) Degree of polynomial f(x)
%
% n_k : (Int) Degree of polynomial v(x)
%
% n : (Int) Degree of polynomial g(x)


% Produce a vector of elements of D, then diagonalise it to form a matrix.
D = diag(1./GetBinomials(m+n_k));

end