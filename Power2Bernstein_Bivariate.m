function [fxy_Bb] = Power2Bernstein_Bivariate(fxy_matrix)


[m1,m2] = GetDegree(fxy_matrix);

A = zeros(m1+1,m1+1);
B = zeros(m2+1,m2+1);

fxy_Bb = A * fxy_matrix * B;

for i = 0:1:m1
    for j = 0:1:m1
        if i >= j
            A(i+1,j+1) = nchoosek(i,j) ./ nchoosek(m1,j)
        end
    end
end

for i = 0:1:m2
    for j = 0:1:m2
        if i >=j
            B(i+1,j+1) = nchoosek(i,j) ./ nchoosek(m2,j);
        end
    end
end

fxy_Bb = A* fxy_matrix * B'



end