function [i,j] = GivenCol_getIndex2(m,n,c)
% m:    Degree with respect to x
% n:    Degree with respect to y
% c:    Index of column removed

% build the matrix of i coefficients
i_matrix = ones(m+1,n+1);

% multiply the rows by 0,1,2,3,4
di_mat = diag(0:1:m);
i_matrix = di_mat * i_matrix ;

i_vec = getAsVector(i_matrix);

% Build the matrix of j coefficients
j_matrix = ones(m+1,n+1);

di_mat = diag(0:1:n);
j_matrix = j_matrix * di_mat;

j_vec = getAsVector(j_matrix);

i = i_vec(c);
j = j_vec(c);