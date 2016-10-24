function G = BuildG(t1,t2)
% Build the matrix G.

G_mat = GetWithBinomials(ones(t1+1,t2+1));

% Get matrix G as a vector
G_vec = GetAsVector(G_mat);

% Diagonalise G
G = diag(G_vec);


end