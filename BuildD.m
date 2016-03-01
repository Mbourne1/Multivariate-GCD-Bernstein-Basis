function D = BuildD(k1,k2,m1,m2,n1,n2)

% initialise the matrix whose entries will form the diagonal entries of D
D_mat = ones(m1+n1-k1+1,m2+n2-k2+1);


% Divide each row of d by the corresponding \binom{m1+n1-k1}{i}
% Divide each col of d by the corresponding \binom{m2+n2-k2}{j}
Bi_rows = GetBinomials(m1+n1-k1);
Bi_cols = GetBinomials(m2+n2-k2);

D_mat = diag(1./Bi_rows) * D_mat * diag(1./Bi_cols);

% Get the D_matrix as a vector
D_vec = GetAsVector(D_mat);

% Form a diagonal matrix from the vector D_vec
D = diag(D_vec);

end

