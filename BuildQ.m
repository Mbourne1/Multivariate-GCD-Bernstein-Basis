
function Q = BuildQ(k1,k2,m1,m2,n1,n2)

Q1 = BuildQ1(n1,n2,k1,k2);
Q2 = BuildQ1(m1,m2,k1,k2);

Q = blkdiag(Q1,Q2);
end