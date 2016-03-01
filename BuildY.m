function Y = BuildY(m1,m2,n1,n2,t1,t2,opt_col,x_ls,alpha,theta1,theta2)

% Note the coefficients of u and v will NOT include binomial coefficients
% if the option bool_q = 'n'.

% insert a zero into the position of the optimal_column
first_part = x_ls(1:(opt_col-1));
second_part = x_ls(opt_col:end);

x = [first_part ; 0 ; second_part];
% partition x_ls into the two parts for coefficients of 


% Separate the x into coefficients of u and coefficients of v

% The first (n1-t1+1) x (n2-t2+1) coefficients are of v
num_coeff_v = (n1-t1+1) * (n2-t2+1);
num_coeff_u = (m1-t1+1) * (m2-t2+1);

% Get the vector of coefficients of v(w,w)
vww_vec = x(1:num_coeff_v);

% Get the vector of coefficients of u(w,w)
uww_vec = - x(num_coeff_v+1:end);

% Get the vector v(w,w) and u(w,w) as a matrix
v_mat = GetAsMatrix(vww_vec,n1-t1,n2-t2);
u_mat = GetAsMatrix(uww_vec,m1-t1,m2-t2);

T2 = BuildC1_new(-u_mat,m1,m2,n1,n2,t1,t2,theta1,theta2);
T1 = BuildC1_new(v_mat,n1,n2,m1,m2,t1,t2,theta1,theta2);


% The last (m1-t1+1) x (m2-t2+1) coefficients are of u

% Multiply by the thetas of f and g
th1_mat = diag(theta1.^(0:1:m1));
th2_mat = diag(theta2.^(0:1:m2));
fww_thetas_mat = th1_mat * ones(m1+1,m2+1) * th2_mat;
fww_thetas_vec = GetAsVector(fww_thetas_mat);

th1_mat = diag(theta1.^(0:1:n1));
th2_mat = diag(theta2.^(0:1:n2));
gww_thetas_mat = th1_mat * ones(n1+1,n2+1) * th2_mat;
gww_thetas_vec = GetAsVector(gww_thetas_mat);

% multiply

vec = [fww_thetas_vec ; gww_thetas_vec];
thetas_mat = diag(vec);


%%
% multiply by the alpha of g
Y = [T1 alpha*T2] * thetas_mat  ;
end