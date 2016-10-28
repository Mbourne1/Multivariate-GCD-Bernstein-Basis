function Y = BuildY(m1,m2,n1,n2,t1,t2,opt_col,x_ls,alpha,th1,th2)
% Build the matrix Y such that D*Y(x)*G*[f;g] = A_{t}(f,g)*x, 
% where A_{t}(f,g) = (D*T(f,g)*Q)*M where M removes an optimal column of
% T(f,g).
%
% Inputs.
%
% m1 :
%
% m2 :
%
% n1 :
%
% n2 :
%
% t1 :
%
% t2 :
%
% opt_col :
%
% x_ls :
%
% alpha :
%
% theta1 :
%
% theta2 :

m1_t1 = m1 - t1;
m2_t2 = m2 - t2;
n1_t1 = n1 - t1; 
n2_t2 = n2 - t2;

% Note the coefficients of u and v will NOT include binomial coefficients
% if the option bool_q = 'n'.

% insert a zero into the position of the optimal_column
first_part = x_ls(1:(opt_col-1));
second_part = x_ls(opt_col:end);

x = [first_part ; 0 ; second_part];
% partition x_ls into the two parts for coefficients of 


% Separate the x into coefficients of u and coefficients of v

% The first (n1-t1+1) x (n2-t2+1) coefficients are of v
num_coeff_v = (n1_t1+1) * (n2_t2+1);


% Get the vector of coefficients of v(w,w)
vww_vec = x(1:num_coeff_v);

% Get the vector of coefficients of u(w,w)
uww_vec = - x(num_coeff_v+1:end);

% Get the vector v(w,w) and u(w,w) as a matrix
vww_mat = GetAsMatrix(vww_vec,n1_t1,n2_t2);
uww_mat = GetAsMatrix(uww_vec,m1_t1,m2_t2);

% Construct T1(v) and T2(u).
T2 = BuildT1(-uww_mat,n1,n2);
T1 = BuildT1(vww_mat,m1,m2);

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