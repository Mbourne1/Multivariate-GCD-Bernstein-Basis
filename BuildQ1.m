function Q1 = BuildQ1(n1,n2,t1,t2)

counter = 1;
for tot = 0:1:(n1-t1)+(n2-t2)
    for i = tot:-1:0
        j = tot-i;
        if i < (n1-t1 +1) && j < (n2-t2 + 1)
            vec(counter) = nchoosek(n1-t1,i) * nchoosek(n2-t2,j);
            counter = counter + 1;
        end
    end
end
Q1 = diag(vec);

end