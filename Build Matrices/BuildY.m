function Y = BuildY(m1,m2,n1,n2,k1,k2,idx_col,x_ls,alpha,th1,th2)
% BuildY(m1,m2,n1,n2,k1,k2,idx_col,x_ls,alpha,th1,th2)
%
%
% Build the matrix Y such that D*Y(x)*G*[f;g] = A_{t}(f,g)*x, 
% where A_{t}(f,g) = (D*T(f,g)*Q)*M where M removes an optimal column of
% T(f,g).
% This function is used in STLN()
%
% Inputs.
%
% m1 : Degree of polynomial f(x,y) with respect to x
%
% m2 : Degree of polynomial f(x,y) with respect to y
%
% n1 : Degree of polynomial g(x,y) with respect to x
%
% n2 : Degree of polynomial g(x,y) with respect to y
%
% k1 : Degree of polynomial d(x,y) with respect to x
%
% k2 : Degree of polynomial d(x,y) with respect to y
%
% idx_col : Index of column removed from Sylvester matrix.
%
% x_ls : 
%
% alpha :
%
% th1 :
%
% th2 :

% Get degree of polynomial u(x,y)
m1_k1 = m1 - k1;
m2_k2 = m2 - k2;

% Get degree of polynomial v(x,y)
n1_k1 = n1 - k1; 
n2_k2 = n2 - k2;


% Insert a zero into the position of the optimal_column.
% Partition x_ls into the two parts for coefficients of 

first_part = x_ls(1:(idx_col-1));
second_part = x_ls(idx_col:end);

x = [first_part ; 0 ; second_part];

% Separate the x into coefficients of u and coefficients of v
% The first (n1-t1+1) x (n2-t2+1) coefficients are of v
nCoeff_v = (n1_k1+1) * (n2_k2+1);


% Get the vector of coefficients of v(w,w)
vww_vec = x(1:nCoeff_v);

% Get the vector of coefficients of u(w,w)
uww_vec = - x(nCoeff_v+1:end);

% Get the vector v(w,w) and u(w,w) as a matrix
vww_mat = GetAsMatrix(vww_vec,n1_k1,n2_k2);
uww_mat = GetAsMatrix(uww_vec,m1_k1,m2_k2);

% Construct T1(v) and T2(u).
T1 = BuildT1(vww_mat,m1,m2);
T2 = BuildT1(-uww_mat,n1,n2);


% Get the matrix of thetas corresponding to (and exlcluding) coefficients
% of f(x,y)
fww_thetas_mat = GetWithThetas(ones(m1+1,m2+1),th1,th2);
fww_thetas_vec = GetAsVector(fww_thetas_mat);

% Get the matrix of thetas corresponding to (and exlcluding) coefficients
% of g(x,y)
gww_thetas_mat = GetWithThetas(ones(n1+1,n2+1),th1,th2);
gww_thetas_vec = GetAsVector(gww_thetas_mat);


vec = [fww_thetas_vec ; gww_thetas_vec];
thetas_mat = diag(vec);

% multiply by the alpha of g
Y = [T1 alpha*T2] * thetas_mat  ;
end