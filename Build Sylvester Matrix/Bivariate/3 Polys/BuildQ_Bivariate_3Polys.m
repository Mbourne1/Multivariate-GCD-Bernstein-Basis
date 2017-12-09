function Q = BuildQ_Bivariate_3Polys(m1, m2, n1, n2, o1, o2, k1, k2)

Q_ok = BuildQ1_Bivariate(o1 - k1, o2 - k2);
Q_nk = BuildQ1_Bivariate(n1 - k1, n2 - k2);
Q_mk = BuildQ1_Bivariate(m1 - k1, m2 - k2);


Q = blkdiag(Q_nk, Q_ok, Q_mk);


end