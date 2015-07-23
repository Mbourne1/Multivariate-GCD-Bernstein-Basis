function exp3 = getInPowerBasis(fxy_matrix_bi)
% Given a matrix fxy_matrix_bi whose entries are
% a_{i,j}\binom{m1}{i}\binom{m2}{j} get an expression in terms of power
% basis

[rows,cols] = size(fxy_matrix_bi);
% Get degree m1 of f (in terms of x)
m1 = rows-1;
% Get degree m2 of f (in terms of y)
m2 = cols-1;

x = sym('x');
y = sym('y');

for i = 0:1:m1
    R(i+1) = x^i * (1-x)^(m1-i);
end

% convert R to a diagonal matrix
R = diag(R);
for j = 0:1:m2
    C(j+1) = y^j * (1-y)^(m2-j);
end

% convert C to a diagonal matrix
C = diag(C);

exp = R * fxy_matrix_bi * C;
exp2 = sum(sum(exp,1),2);
exp3 = expand(exp2);

end