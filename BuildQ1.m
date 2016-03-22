function Q1 = BuildQ1(n1_t1,n2_t2)
% Build the matrix Q1, a partition of the matrix Q, where Q is the diagonal
% matrix of binomial coefficients of the vector [f;g].

Q_mat = GetWithBinomials(ones(n1_t1+1,n2_t2+1));

% Get the matrix Q as a vector
Q_vec = GetAsVector(Q_mat);

% Diagonalise Q
Q1 = diag(Q_vec);

end