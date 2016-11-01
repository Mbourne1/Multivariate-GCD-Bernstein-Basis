function D = BuildD(m1,m2,n1_k1,n2_k2)
% Build the matrix D

% Divide each row of d by the corresponding \binom{m1+n1-k1}{i}
% Divide each col of d by the corresponding \binom{m2+n2-k2}{j}
D_mat = GetWithBinomials(ones(m1+n1_k1+1,m2+n2_k2+1));

% Get the D_matrix as a vector
D_vec = GetAsVector(1./D_mat);

% Form a diagonal matrix from the vector D_vec
D = diag(D_vec);

end

