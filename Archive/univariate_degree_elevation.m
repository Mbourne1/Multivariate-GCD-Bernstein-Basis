function [f_new] = univariate_degree_elevation(fx)
% fx = coefficients of fx
% r = number to increase degree by

m = length(fx)-1;

f_new = zeros(1,m+2);

for i = 0:1:m+1
    if i == 0
        f_new(i+1) = fx(i+1) * (m+1-i)./(m+1);
    elseif i == m+1
        f_new(i+1) = fx(i) * (i)./(m+1);
    else
        f_new(i+1) =...
            fx(i) * (i)./(m+1) ...
            +...
            fx(i+1) * (m+1-i)./(m+1);
    end
end

r = 1;
c = zeros(1,m+r+1);

for k = 0:1:m+r
    c(k+1) = summing_function(fx,r,k,m);
end

f_new
c

end

function sum = summing_function(fx,r,k,m)
sum = 0;
for j = max(0,k-r):1:min(m,k)
    sum = sum + (nchoosek(r,k-j) * nchoosek(m,j) ./ nchoosek(m+r,k) * fx(j+1));
end
end
