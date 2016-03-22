function Q = BuildQ(n1_k1,n2_k2,m1_k1,m2_k2)
% Build the diagonal matrix Q whose entries contain binomial coefficients
% corresponding to the entries of u(x,y) and v(x,y) in the vector [v;-u]

Q1 = BuildQ1(n1_k1,n2_k2);
Q2 = BuildQ1(m1_k1,m2_k2);

Q = blkdiag(Q1,Q2);

end