function Q1 = BuildQ1(n1,n2,t1,t2)

% Build the matrix Q
Q_mat = ones(n1-t1+1,n2-t2+1);

% Multiply rows by binomials
% Multiply columns by binomials
Bi_n1_t1 = GetBinomials(n1-t1);
Bi_n2_t2 = GetBinomials(n2-t2);

Q_mat = diag(Bi_n1_t1) * Q_mat * diag(Bi_n2_t2);

% Get the matrix Q as a vector
Q_vec = GetAsVector(Q_mat);

% Diagonalise Q
Q1 = diag(Q_vec);




end