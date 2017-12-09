function Q1 = BuildQ1_Bivariate(n1_k1, n2_k2)
% BuildQ1_Bivariate(n1_t1, n2_t2)
%
% Build the matrix Q1, a partition of the matrix Q, where Q is the diagonal
% matrix of binomial coefficients of the vector [v;-u].
% Q is a matrix which forms the kth Sylvester subresultant matrix
% D^{-1}[T(f) T(g)]Q
%
% Inputs
%
% n1_k1 : (Int) Degree of v(x,y) with respect to x
%
% n2_k2 : (Int) Degree of v(x,y) with respect to y
%
% Outputs.
%
% Q1 : (Matrix) Partition Q_{1} of the matrix Q. Such that Q_{1} * v where v is a
% vector of coefficients of v(x,y) gives the coefficients in scaled
% Bernstein form.

% Get a matrix of binomial coefficients corresponding to the polynomial
% v(x,y).
Q_mat = GetWithBinomials_Bivariate(ones(n1_k1 + 1, n2_k2 + 1));


% Get the matrix Q as a vector
Q_vec = GetAsVector_Version1(Q_mat);


% Diagonalise Q
Q1 = diag(Q_vec);

end