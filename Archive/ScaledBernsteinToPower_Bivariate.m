function fxy = ScaledBernsteinToPower_Bivariate(fxy_matrix_bi)
% Given a matrix fxy_matrix_bi whose entries are
% fxy_matrix_bi(i+1,j+1) = a_{i,j} \binom{m1}{i}\binom{m2}{j} 
% get the polynomial in the power basis.

% Get the degree of fxy wrt x and y
[m1, m2] = GetDegree(fxy_matrix_bi);

% %
% Initialise some symbolic variables
x = sym('x');
y = sym('y');


% % Build Matrix R

% for each row
R = zeros(m+1,1);
for i = 0:1:m1
    R(i+1) = x^i * (1-x)^(m1-i);
end

% Convert R to a diagonal matrix
R = diag(R);

% % Build Matrix C

% for each column
for j = 0:1:m2
    C(j+1) = y^j * (1-y)^(m2-j);
end

% Convert C to a diagonal matrix
C = diag(C);

% % Multiply R * f * C to obtain symbolic expression of f in power basis

exp = R * fxy_matrix_bi * C;

exp2 = sum(sum(exp,1),2);

fxy = expand(exp2);

end