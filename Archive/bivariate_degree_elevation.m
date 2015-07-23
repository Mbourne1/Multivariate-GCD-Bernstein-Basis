function [ F ] = bivariate_degree_elevation( fxy_matrix,r,s )
%BIVARIATE_DEGREE_ELEVATION Summary of this function goes here
%   Detailed explanation goes here
% r = degree elevation corresponding to m
% s = degree elevation corresponding to n
%

% given an input matrix where the rows are in the second basis (y) and
% columns are first basis (x)
% The entries exclude binomial coefficients

fxy_matrix = fxy_matrix';

[m1,n1] = size(fxy_matrix);
n = n1-1;
m = m1-1;


% for each column j
for j = 0:1:n+s
    % for each row i
    for i = 0:1:m+r
        F(i+1,j+1) = summing_function(fxy_matrix,m,n,r,s,i,j);
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