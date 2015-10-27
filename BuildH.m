function H = BuildH(m1,m2,n1,n2)


    H1 = BuildH1(m1,m2);
    H2 = BuildH1(n1,n2);

    H = blkdiag(H1,H2);
end