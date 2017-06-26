function P = BuildP1(m1,m2,n1,n2,k1,k2,idx_col)
% BuildP1(m1,m2,n1,n2,th1,th2,idx_col,k1,k2)
%
% Build the matrix P1 which forms a partition of the matrix P, and is used
% in the SNTLN function. See BuildP() for more details
%
% Inputs.
%
% m1 : (Int) Degree of polynomial f(x,y) with respect to x
%
% m2 : (Int) Degree of polynomial f(x,y) with respect to y
%
% n1 : (Int) Degree of polynomial g(x,y) with respect to x
%
% n2 : (Int) Degree of polynomial g(x,y) with respect to y
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}
%
% k1 : (Int) Degree of polynomial d(x,y) with respect to x
%
% k2 : (Int) Degree of polynomial d(x,y) with respect to y
%
% idx_col : (Int) Index of column removed from the Sylvester Subresultant
% matrix S_{k_{1},k_{2}}(f,g) such that A_{k_{1},k_{2}}x = c_{k_{1},k_{2}}.


% Get a matrix of \theta_{1}^{i}\theta_{2}^{j} corresponding to polynomial
% f(x,y). 
thetas_matrix = ones(m1+1, m2+1);

% Produce a zero matrix to fill the space
padd_mat = zeros(m1 + n1 - k1 + 1, m2 + n2 - k2 + 1);

% % Given the index of the optimal column, Get the number of multiplications
% with respec to x and number with respect to y.
[i,j] = GivenCol_GetIndex(n1 - k1, n2 - k2, idx_col);

% initialise i hat and j hat.
ihat = i+1;
jhat = j+1;

% Get the number of rows in f(x,y)
nRows_f = m1 + 1;

% Get the number of columns in f(x,y)
nColumns_f = m2 + 1;

% Place the thetas matrix (matrix of thetas corresponding to f(x,y)) into
% the padded matrix.
padd_mat(ihat : i + nRows_f, jhat : j + nColumns_f) = thetas_matrix;

% Get the matrix as a vector, starting in upper left corner, read down and
% diagonally upwards.
vec_padd_mat = GetAsVector_Version1(padd_mat);

% Form a diagonal matrix of the vector we just produced.
diag_mat_vec_padd_mat = diag(vec_padd_mat);

% Remove the zero columns
diag_mat_vec_padd_mat( :, ~any(diag_mat_vec_padd_mat,1) ) = [];  %columns

% Assign to P
P = diag_mat_vec_padd_mat;

P = nchoosek(n1 - k1, i) * nchoosek(n2 - k2, j) * P;

end