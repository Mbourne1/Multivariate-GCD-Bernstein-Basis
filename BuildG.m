function G = BuildG(t1,t2)


G_mat = ones(t1+1,t2+1);

% Multiply each row of G by corresponding binomial coefficient
for i=0:1:t1
   G_mat(i+1,:) =  G_mat(i+1,:) .* nchoosek(t1,i);
end

% Multiply each column of G by corresponding binomial coefficient
for j=0:1:t2
    G_mat(:,j+1) = G_mat(:,j+1) .* nchoosek(t2,j);
end

G_vec = getAsVector(G_mat);

G = diag(G_vec);


end