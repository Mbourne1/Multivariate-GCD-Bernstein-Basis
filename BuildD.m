function D = BuildD(k1,k2,m1,m2,n1,n2)

% initialise the matrix whose entries will form the diagonal entries of D
D_mat = ones(m1+n1-k1+1,m2+n2-k2+1);

% Divide each row of d by the corresponding \binom{m1+n1-k1}{i}
for i = 0:1:m1+n1-k1
   D_mat(i+1,:) = D_mat(i+1,:) ./ nchoosek(m1+n1-k1,i); 
end

% Divide each col of d by the corresponding \binom{m2+n2-k2}{j}
for j = 0:1:m2+n2-k2
   D_mat(:,j+1) = D_mat(:,j+1) ./ nchoosek(m2+n2-k2,j); 
end

% Get the D_matrix as a vector
D_vec = getAsVector(D_mat);

% Form a diagonal matrix from the vector D_vec
D = diag(D_vec);

end

