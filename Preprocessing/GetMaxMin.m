function [max_matrix, min_matrix] = GetMaxMin(fxy, n1_k1, n2_k2)
% GetMaxMin(fxy, n1_k1, n2_k2)
%
% Get the maximum and minimum of each entry of f(x,y) in the Sylvester
% matrix.
%
% Inputs
%
% fxy : (Matrix) The Coefficients of polynomial f(x,y) in standard bernstein basis.
%       Given in matrix form so that the rows are in terms of x basis
%       elements and the columns are y basis elements.
%
% n1_k1 : (Int) Degree of v(x,y) with respect to x
%
% n2_k2 : (Int) Degree of v(x,y) with respect to y
%
% % Outputs
%
% max_matrix : (Matrix)
%
% min_matrix : (Matrix)

% Get the degree of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Take the absolute values of the coefficients
fxy = abs(fxy);

% Build a matrix which stores the maximum and minimum values of each coefficient
max_matrix = zeros(m1 + 1, m2 + 1);
min_matrix = zeros(m1 + 1, m2 + 1);

for i1=0 : 1 : m1
    for i2 = 0 : 1 : m2
        
        % This gives the a_{i1,i2} coefficient of f
        %[maximum,minimum] = GetMaxMin2(fxy(i1+1,i2+1),i1,i2,m1,m2,n1_k1,n2_k2);

        [maximum, minimum] = GetMaxMin2a(fxy(i1+1, i2+1), i1, i2, m1, m2, n1_k1, n2_k2);
        max_matrix(i1 + 1, i2 + 1) = maximum;
        min_matrix(i1 + 1, i2 + 1) = minimum;
        
    end
    
end


end


% % The next three functions compute max and minimum entries by three
% different methods.

function [maximum,minimum] = GetMaxMin2a(a_i1i2, i1, i2, m1, m2, n1_k1, n2_k2)
% 
%
%

% Build a 2 dimensional vector to store all occurences of the coefficient
% a_{i_{1},i_{2}}

val = abs(a_i1i2) .* nchoosek(m1, i1) * nchoosek(m2, i2);

A = ones(n1_k1 + 1, n2_k2 + 1) .* val;

B1 = zeros(n1_k1 + 1, 1);
B2 = zeros(n1_k1 + 1, 1);

C1 = zeros(n2_k2 + 1, 1);
C2 = zeros(n2_k2 + 1, 1);


for j1 = 0 : 1 : n1_k1
    B1(j1 + 1) = nchoosek(n1_k1, j1);
    B2(j1 + 1) = nchoosek(m1 + n1_k1, i1 + j1);
end

for j2 = 0 : 1 : n2_k2
    C1(j2+1) = nchoosek(n2_k2, j2);
    C2(j2+1) = nchoosek(m2 + n2_k2, i2 + j2);
end

global SETTINGS

switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    
    case 'T'
        matrix = A;
        
    case 'DT'
        matrix = diag(1./B2) * A * diag(1./C2);
        
    case 'TQ'
        matrix = diag(B1) * A * diag(C1);
        
    case 'DTQ'
        matrix = diag(B1./B2) * A * diag(C1./C2);
        
    case 'DTQ Denominator Removed'
        matrix = diag(B1./B2) * A * diag(C1./C2) .* nchoosek(m1 + n1_k1, n1_k1) .* nchoosek(m2 + n2_k2, n2_k2);
        
    otherwise 
        error('err')
end






% take absolute values of A
A = abs(matrix);

[max_row, max_col] = find(A==max(A(:)));

[min_row, min_col] = find(A==min(A(:)));


% get the maximum and minimum values. Always use (1) since max or min may
% occur more than once, and we are only interested in one of these values.
maximum = A(max_row(1),max_col(1));
minimum = A(min_row(1),min_col(1));



end
