function G = BuildG(t1,t2)
% Build the matrix G.

% Initialise matrix G
G_mat = ones(t1+1,t2+1);

% Multiply rows by nchoosek(t1,i)
% Multiply columns by nchoosek(t2,j)
Bi_t1 = GetBinomials(t1);
Bi_t2 = GetBinomials(t2);
G_mat = diag(Bi_t1) * G_mat * diag(Bi_t2);

% Get matrix G as a vector
G_vec = GetAsVector(G_mat);

% Diagonalise G
G = diag(G_vec);


end