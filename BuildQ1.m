function Q1 = BuildQ1(n1,n2,t1,t2)

% Build the matrix Q

Q_mat = ones(n1-t1+1,n2-t2+1);

% Multiply rows by binomials
for i = 0:1:n1-t1
   Q_mat(i+1,:) = Q_mat(i+1,:) .* nchoosek(n1-t1,i) ;
end

% Multiply cols by binomials
for j = 0:1:n2-t2
   Q_mat(:,j+1) = Q_mat(:,j+1) .* nchoosek(n2-t2,j) ;
end

Q_vec = getAsVector(Q_mat);
Q1 = diag(Q_vec);




end