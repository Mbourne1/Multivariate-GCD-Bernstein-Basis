% This file will test the scaling of the binomial coefficients


m1 = 5;
m2 = 4;

n1_k1 = 3;
n2_k2 = 4;


fxy = ones(m1+1,m2+1);

mySum = zeros(m1+1, m2+1);

for i1 = 0:1:m1
    for i2 = 0:1:m2
        
        coefficient_sum = 0;

        for j1 = 0:1:n1_k1
            for j2 = 0:1:n2_k2

                temp_val = ...
                    nchoosek(i1 + j2, i1) ...
                    * nchoosek(i2+j2, i2) ...
                    * nchoosek(m1 + n1_k1 - (i1+j1), m1-i1) ...
                    * nchoosek(m2 + n2_k2 - (i2+j2), m2-i2) ...
                    / (nchoosek(m1 + n1_k1, m1) * nchoosek(m2 + n2_k2, m2) );
                
                coefficient_sum = coefficient_sum + temp_val;

            end
        end

        sum_matrix(i1+1,i2+1) = coefficient_sum;

    end
end




