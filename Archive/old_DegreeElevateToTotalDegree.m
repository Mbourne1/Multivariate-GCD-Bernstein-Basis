function deg_elv_fxy_matrix = old_DegreeElevateToTotalDegree(fxy_matrix,m)


[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% let r_f be the number of rows being added
r_f = m - m1; 

% let s_f be the number of cols being added
s_f = m - m2;

% Degree elevate fxy_matrix
deg_elv_fxy_matrix = zeros(m+1,m+1);

% for each column j
for j = 0:1:m2 + s_f
    % for each row i
    for i = 0:1:m1 + r_f
        deg_elv_fxy_matrix(i+1,j+1) = summing_function(fxy_matrix,m1,m2,r_f,s_f,i,j);
    end
end


end

function sum = summing_function(fxy_matrix,m,n,r,s,i,j)

sum = 0;
for k = max(0,i-r) : 1 : min(m,i)
    for l = max(0,j-s) : 1 : min(n,j)
        sum = sum +...
            (...
            fxy_matrix(k+1,l+1) *...
            nchoosek(m,k) *...
            nchoosek(r,i-k) *...
            nchoosek(n,l) *...
            nchoosek(s,j-l) /...
            ...
            (nchoosek(m+r,i) *...
            nchoosek(n+s,j))...
            );
    end
end

end

