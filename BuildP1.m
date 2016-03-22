function P = BuildP1(m1,m2,n1,n2,th1,th2,opt_col,t1,t2)
% Build the matrix P such that P*[f;g] = c_{t}, where c_{t} is an optimal
% column removed from the Sylvester matrix.

% Get a matrix of thetas corresponding to polynomial f(x,y)
thetas_matrix = GetWithThetas(ones(m1+1,m2+1),th1,th2)

% Produce a zero matrix to fill the space
padd_mat = zeros(m1+n1-t1+1, m2+n2-t2+1);

% % Given the index of the optimal column, Get the number of multiplications
% with respec to x and number with respect to y.
[i,j] = GivenCol_getIndex2(n1-t1,n2-t2,opt_col);

% initialise i hat and j hat.
ihat = i+1;
jhat = j+1;

% Get the number of rows in f(x,y)
num_rows_f = m1+1;

% Get the number of columns in f(x,y)
num_cols_f = m2+1;

% Place the thetas matrix (matrix of thetas corresponding to f(x,y)) into
% the padded matrix.
padd_mat(ihat:i+num_rows_f, jhat:j+num_cols_f) = thetas_matrix;

% Get the matrix as a vector, starting in upper left corner, read down and
% diagonally upwards.
vec_padd_mat = GetAsVector(padd_mat);

% Form a diagonal matrix of the vector we just produced.
diag_mat_vec_padd_mat = diag(vec_padd_mat);


% Remove the zero columns
diag_mat_vec_padd_mat( :, ~any(diag_mat_vec_padd_mat,1) ) = [];  %columns

% Assign to P
P = diag_mat_vec_padd_mat;

P = nchoosek(n1-t1,i) * nchoosek(n2-t2,j) * P;

end