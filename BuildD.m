function D = BuildD(k1,k2,m1,m2,n1,n2)

count =1 ;

for tot = 0:1:m1+m2+n1+n2-k1-k2
    for i_hat = tot:-1:0
        j_hat = tot - i_hat;
        if j_hat < m2+n2-k2+1 && i_hat < m1+n1-k1+1
            [i_hat,j_hat];
            temp_vec(count) = nchoosek(m1+n1-k1,i_hat) * nchoosek(m2+n2-k2,j_hat);
            count = count+1;
        end
    end
end
D = diag(1./temp_vec);
end

