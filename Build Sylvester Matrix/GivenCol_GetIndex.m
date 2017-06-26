function [i,j] = GivenCol_GetIndex(m,n,c)
% GivenCol_GetIndex(m,n,c)
%
%
% % Inputs.
%
% m : (Int) Degree with respect to x
%
% n : (Int) Degree with respect to y
%
% c : (int)   Index of column removed
%
% % Outputs
%
% [i,j] : Index of column

% build the matrix of i coefficients
i_matrix = ones(m+1, n+1);

% multiply the rows by 0,1,2,3,4
di_mat = diag(0 : 1 : m);
i_matrix = di_mat * i_matrix ;

i_vec = GetAsVector_Version1(i_matrix);

% Build the matrix of j coefficients
j_matrix = ones(m+1, n+1);

di_mat = diag(0 : 1 : n);
j_matrix = j_matrix * di_mat;

j_vec = GetAsVector_Version1(j_matrix);

i = i_vec(c);
j = j_vec(c);

end