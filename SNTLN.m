function [ fxy_output,gxy_output,alpha_output,theta1_output,theta2_output,X_output] = ...
    SNTLN( fxy_matrix,gxy_matrix,...
    initial_alpha, initial_theta_1, initial_theta_2,...
    t1,t2,...
    lambda,mu,...
    opt_col)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                             Inputs.


% fx_n :    Coefficients of polynomial f, in standard bernstein basis.

% gx_n :    Coefficients of polynomial g, in standard bernstein basis.

% initial_alpha :   Initial value of alpha

% initial_theta :   Initial value of theta

% t :   Degree of AGCD.

% opt_col : Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Outputs.


% fx_output :- Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.

% gx_output :- Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.

% alpha_output :-

% theta_output :-

% X_output :-

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Inputs

global max_error
global max_iterations

max_error = 1e-14;
max_iterations = 100;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the initial iterations number
ite = 1;

% Set initial values of alpha and theta
theta1(ite) = initial_theta_1;
theta2(ite) = initial_theta_2;
alpha(ite) = initial_alpha;

% Get degree of polynomials f.
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get degree of polynomial g.
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

% Get the number of coefficients in the polynomial f(x,y)
num_coeff_f = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
num_coeff_g = (n1+1) * (n2+1);

% Get the number of coefficients in both f(x,y) and g(x,y)
num_coeff = num_coeff_f + num_coeff_g;

% Get the number of coefficients in v(x,y)
num_coeff_v = (n1-t1+1) * (n2-t2+1);

% Get the number of coefficients in u(x,y)
num_coeff_u = (m1-t1+1) * (m2-t2+1);

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
num_coeff_x = num_coeff_u + num_coeff_v - 1;

% Create the identity matrix I, the matrix M formed from I by removing the
% column equivalent to the optimal column for removal from the Sylvester
% subresultant matrix, so that S_{t}(f,g)*M = A_{t}, where A_{t} is the
% Sylvester subresultant matrix with the column removed.

% Get the number of columns in C_{t}(f), the first partition of the Sylvester
% Matrix S_{t}(f,g)
num_cols_in_f_partition = (n1-t1+1) * (n2-t2+1);

% Get the number of columns in C_{t}(g), the second partition of the
% Sylvester matrix S_{t}(f,g)
num_cols_in_g_partition = (m1-t1+1) * (m2-t2+1);

% Get the total number of columns in the Sylvester matrix S_{t}(f,g)
num_cols = num_cols_in_f_partition + num_cols_in_g_partition;

% Create the identity matrix
I = eye(num_cols, num_cols);

% Create the matrix M, such that S(f,g)*M gives A_{t}, the Sylvester Matrix
% with the optimal column removed.
M = I;
M(:,opt_col) = [];

% Let e be the column removed from the identity matrix, such that
% S_{t}(f,g) * e gives the column c_{t}, where c_{t} is the optimal column
% removed from the Sylvester subresultant.
e = I(:,opt_col);

%% Preprocessing

% Normalise polynomial f(x,y) by geometric means
fxy_matrix_n = fxy_matrix./lambda;
gxy_matrix_n = gxy_matrix./mu;

% Obtain polynomials in Modified Bernstein Basis, using initial values of
% alpha and theta.

% Multiply the rows of fxy_matrix by theta1, and multiply the cols of
% fxy_matrix by theta2
fw1w2_matrix = diag(theta1.^(0:1:m1)) * fxy_matrix_n * diag(theta2.^(0:1:m2));

% Multiply the rows of gxy_matrix by theta1, and multiply the cols of
% gxy_matrix by theta2
gw1w2_matrix = diag(theta1.^(0:1:n1)) * gxy_matrix_n * diag(theta2.^(0:1:n2));

%% Form the Coefficient Matrix T = [C(f)|C(g)] such that DTQ * x = [col]

D = BuildD(t1,t2,m1,m2,n1,n2);
Q = BuildQ(t1,t2,m1,m2,n1,n2);

T1 = BuildT1(fw1w2_matrix,n1,n2,t1,t2,1,1);
T2 = BuildT1(gw1w2_matrix,m1,m2,t1,t2,1,1);
DTQ = D*[T1 alpha(ite).*T2]*Q;


%%
% Calculate the partial derivatives of f(w,w) and g(w,w) with respect to alpha
Partial_fw_wrt_alpha            = zeros(m1+1,m2+1);
Partial_alpha_gw_wrt_alpha      = gxy_matrix_n;

%%
% Calculate the partial derivatives of f(w,w) with respect to theta_1
theta_mat = diag((0:1:m1)./theta1(ite));
Partial_fw_wrt_theta1    = theta_mat * fw1w2_matrix;

%%
% Calculate the partial derivative of g(w,w) with respect to theta_1
theta_mat = diag((0:1:n1)./theta1(ite));
Partial_gw_wrt_theta1    = theta_mat * gw1w2_matrix;

%%
% Calculate the partial derivative of f(w,w) with respect to theta_2
theta_mat = diag((0:1:m2)./theta2(ite));
Partial_fw_wrt_theta2 = fw1w2_matrix * theta_mat;

%%
% Calculate the partial deriviates of g(w,w) with respect to theta_2
theta_mat = diag((0:1:n2)./theta2(ite));
Partial_gw_wrt_theta2 = gw1w2_matrix * theta_mat;

%%
% Build the derivative of T(f,g) with respect to alpha
T1_wrt_alpha = BuildT1(Partial_fw_wrt_alpha,n1,n2,t1,t2,1,1);
T2_wrt_alpha = BuildT1(Partial_alpha_gw_wrt_alpha,m1,m2,t1,t2,1,1);
DTQ_alpha = D*[T1_wrt_alpha T2_wrt_alpha]*Q;

%%
% Calculate the derivative of T(f,g) with respect to theta_{1}
T1_wrt_theta1 = BuildT1(Partial_fw_wrt_theta1,n1,n2,t1,t2,1,1);
T2_wrt_theta1 = BuildT1(Partial_gw_wrt_theta1,m1,m2,t1,t2,1,1);
Partial_DTQ_wrt_theta1 = D*[T1_wrt_theta1 alpha(ite)* T2_wrt_theta1]*Q;

%%
% Calcualte the derivative of T(f,g) with respect to theta_2
T1_wrt_theta2 = BuildT1(Partial_fw_wrt_theta2,n1,n2,t1,t2,1,1);
T2_wrt_theta2 = BuildT1(Partial_gw_wrt_theta2,m1,m2,t1,t2,1,1);
Partial_DTQ_wrt_theta2 = D*[T1_wrt_theta2 alpha(ite)*T2_wrt_theta2]*Q;

%%
% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.

num_rows_Sylv_mat = (m1 + n1 - t1 + 1) * (m2 + n2 - t2 + 1);
num_cols_T1 = (n1 - t1 + 1) * (n2 - t2 + 1);
num_cols_T2 = (m1 - t1 + 1) * (m2 - t2 + 1);

zk = zeros(num_coeff , 1);

%%
% Initilaise the derivative of N wrt alpha.
Partial_DNQ_wrt_alpha   = zeros(num_rows_Sylv_mat,num_cols_T1 + num_cols_T2);

% Initilaise the derivative of N wrt theta_1.
Partial_DNQ_wrt_theta1   = zeros(num_rows_Sylv_mat,num_cols_T1 + num_cols_T2);

% Initialise the derivative of N wrt theta 2
Partial_DNQ_wrt_theta2   = zeros(num_rows_Sylv_mat,num_cols_T1 + num_cols_T2);

%%
% Initialise the derivative of h 
% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
Partial_h_wrt_alpha     = Partial_DNQ_wrt_alpha*e;
Partial_h_wrt_theta1    = Partial_DNQ_wrt_theta1*e;
Partial_h_wrt_theta2    = Partial_DNQ_wrt_theta2*e;

% Get the coefficients of f(x,y) in matrix form
fxy_vec = getAsVector(fxy_matrix_n);
gxy_vec = getAsVector(gxy_matrix_n);

% Get the matrix A_{k}(f,g), which is the subresultant matrix S(f,g) with
% an opitmal column removed
Ak = DTQ;
ck = DTQ(:,opt_col);
Ak(:,opt_col) = [];

%%
% % Build the matrix of binomials corresponding to polynomial f(x,y)
binoms_f = ones(m1+1,m2+1);
% Multiply each row by nchoosek(m1,i)
for i = 0:1:m1
   binoms_f(i+1,:) = binoms_f(i+1,:) * nchoosek(m1,i);
end
% Multiply each col by nchoosek(m2,j)
for j = 0:1:m2
    binoms_f(:,j+1) = binoms_f(:,j+1) * nchoosek(m2,j);
end

% % Build the matrix of binomials corresponding to polynomial g(x,y)
binoms_g = ones(n1+1,n2+1);
% Multiply each row by nchoosek(n1,i)
for i = 0:1:n1
   binoms_g(i+1,:) = binoms_g(i+1,:) * nchoosek(n1,i);
end
% Multiply each col by nchoosek(m2,j)
for j = 0:1:n2
    binoms_g(:,j+1) = binoms_g(:,j+1) * nchoosek(n2,j);
end

binoms_f_vec = getAsVector(binoms_f);
binoms_g_vec = getAsVector(binoms_g);
G = diag([binoms_f_vec;binoms_g_vec]);
%%

P = BuildP(m1,m2,n1,n2,theta1(ite),theta2(ite),alpha(ite),t1,t2,opt_col);
P = D*P*G;

%%
% Calculate the derivatives wrt alpha and theta of the removed column.
Partial_ck_wrt_alpha        = DTQ_alpha*e;
Partial_ck_wrt_theta1       = Partial_DTQ_wrt_theta1*e;
Partial_ck_wrt_theta2       = Partial_DTQ_wrt_theta2*e;


%%
% Perform QR decomposition of Ak to obtain the solution x
[~,n_col] = size(Ak);
[Qr,R] = qr(Ak);
R1 = R(1:n_col,:);
cd = Qr'*ck;
c = cd(1:n_col,:);
x_ls = R1\c;

first_part = x_ls(1:(opt_col-1));
second_part = x_ls(opt_col:end);
x = [first_part ; 0 ; second_part];

% Build Matrix Y, where Y(v,u)*[f;g] = S(f,g)*[u;v]
Y = BuildY(m1,m2,n1,n2,t1,t2,opt_col,x_ls,alpha(ite),theta1(ite),theta2(ite));
DYG = D*Y*G;

test1 = DYG * [fxy_vec;gxy_vec];
test2 = DTQ * x;
test1./test2;

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ck - (DTQ*M*x_ls);

%% Get the matrix p, which will store all the perturbations returned from LSE file
num_entries = num_coeff_f...
    + num_coeff_g ...
    + num_coeff_x ...
    + 3;


% Set the initial value of vector p to be zero
p = zeros(num_coeff_f...
    + num_coeff_g ...
    + num_coeff_x ...
    + 3,1);
%%
% Set the intial value of E to the identity matrix
E = eye(num_entries);

%%
% Create the matrix D(T+N)Q, initially N is empty so this is the same as T.
%T1 = BuildT1(fxy_matrix_n,n1,n2,t1,t2,1,1);
%T2 = BuildT1(gxy_matrix_n,m1,m2,t1,t2,1,1);

% Edit 17/11/2015 15:00
T1 = BuildT1(fw1w2_matrix,n1,n2,t1,t2,1,1);
T2 = BuildT1(gw1w2_matrix,m1,m2,t1,t2,1,1);

DTNQ = D*[T1 alpha(ite).*T2]*Q;

%%
% Create The matrix (T+N) with respect to alpha
T1 = BuildT1(Partial_fw_wrt_alpha,n1,n2,t1,t2,1,1);
T2 = BuildT1(Partial_alpha_gw_wrt_alpha,m1,m2,t1,t2,1,1);
DTNQ_wrt_alpha = D*[T1 T2]*Q;

%%
% Create The matrix (T+N) with respect to theta1
T1 = BuildT1(Partial_fw_wrt_theta1,n1,n2,t1,t2,1,1);
T2 = BuildT1(Partial_gw_wrt_theta1,m1,m2,t1,t2,1,1);
DTNQ_wrt_theta1 = D*[T1 alpha(ite) * T2]*Q;

%%
% Create The matrix (T+N) with respect to theta2
T1 = BuildT1(Partial_fw_wrt_theta2,n1,n2,t1,t2,1,1);
T2 = BuildT1(Partial_gw_wrt_theta2,m1,m2,t1,t2,1,1);
DTNQ_wrt_theta2 = D*[T1 alpha(ite) * T2]*Q;

%%
% Create the matrix C for input into iteration

H_z     = DYG-P;

H_x     = DTNQ*M;

H_alpha  = DTNQ_wrt_alpha*M*x_ls - ...
    (Partial_ck_wrt_alpha + Partial_h_wrt_alpha);

H_theta1 = DTNQ_wrt_theta1*M*x_ls - ...
    (Partial_ck_wrt_theta1 + Partial_h_wrt_theta1);

H_theta2 = DTNQ_wrt_theta2*M*x_ls - ...
    (Partial_ck_wrt_theta2 + Partial_h_wrt_theta2);

C       = [H_z H_x H_alpha H_theta1 H_theta2];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    zk;...
    x_ls;...
    alpha(ite);...
    theta1(ite);...
    theta2(ite)
    ];

yy              =   start_point;

% Set the termination criterion to a large value. It will be
% over written later.

%condition(ite) = norm(res_vec)/norm(ck);
% Edit 17/11/2015
condition(ite) = norm(res_vec);
xk = x_ls;

%%%%%%%%%%%%%%%%%%%%%%%%%%%


while condition(ite) >(max_error) &&  ite < max_iterations
    %while   ite < max_iterations
    
    
    % Use the QR decomposition to solve the LSE problem
    % min |y-p| subject to Cy=q
    
    y = LSE(E,p,C,res_vec);
    y_old = y;
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % obtain the small changes

    %% Break down y into its sections
    
    % get the coefficients corresponding to f and g
    delta_zk        = y(1:num_coeff_f + num_coeff_g ,1);
    
    % Remove the zk coefficients from the list of coefficients
    y(1:num_coeff_f + num_coeff_g) = [];
    
    % Get the coefficients corresponding to x
    delta_xk        = y(1:num_coeff_x,1);
    
    % Remove them from the list of coefficients
    y(1:num_coeff_x) = [];
    
    % Get the coefficient corresponding to alpha
    delta_alpha     = y(1:1);
    y(1) = [];
    
    % Get the coefficient corresponding to theta1
    delta_theta1    = y(1:1);
    y(1) = [];
    
    % Get the coefficient corresponding to theta2
    delta_theta2    = y(1:1);
    y(1) = [];
    
    %% Update the variables
    
    % Update variables z_{k}, where z_{k} are perturbations in the
    % coefficients of f and g. 
    zk = zk + delta_zk;
    
    % Update x_{k}, where x_{k} is the solution vector, containing
    % coefficients u and v.
    xk = xk + delta_xk;
    
    % Update alpha 
    alpha(ite) = alpha(ite-1) + delta_alpha;
    
    % Update theta_{1}
    theta1(ite) = theta1(ite-1) + delta_theta1;
    
    % Update theta_{2}
    theta2(ite) = theta2(ite-1) + delta_theta2;
    
    % Obtain polynomials in modified bersntein basis a_{i}\theta^{i}
    
    % Obtain new f(w,w) with improved theta1, and theta2
    pre_theta1      = diag(theta1(ite).^(0:1:m1));
    post_theta2     = diag(theta2(ite).^(0:1:m2));
    fw1w2_matrix    = pre_theta1 * fxy_matrix_n * post_theta2;
    
    % Obtain new g(w,w) with improved theta1 and theta2
    pre_theta1      = diag(theta1(ite).^(0:1:n1));
    post_theta2     = diag(theta2(ite).^(0:1:n2));
    gw1w2_matrix    = pre_theta1 * gxy_matrix_n * post_theta2;
    
    
    % Construct the Sylvester subresultant matrix DTQ.   
    T1 = BuildT1(fw1w2_matrix,n1,n2,t1,t2,1,1);
    T2 = BuildT1(gw1w2_matrix,m1,m2,t1,t2,1,1);
    DTQ = D*[T1 alpha(ite).*T2]*Q;
    
    % Calculate the partial derivatives of fw and gw with respect to alpha
    Partial_fw_wrt_alpha            = zeros(m1+1,m2+1);
    Partial_alpha_gw_wrt_alpha      = gw1w2_matrix;
    
    %%
    % Calculate the partial derivatives of fw and gw with respect to theta1
    % divide the rows by theta1 and multiply by the old power
    
    % Get the partial derivative of f(w1,w2) with respect to theta_{1}
    temp_mat = diag((0:1:m1)./theta1(ite));
    Partial_fw_wrt_theta1 = temp_mat * fw1w2_matrix;
    
    % Get the partial derivative of g(w1,w2) with respect to theta_{1}
    temp_mat = diag((0:1:n1)./theta1(ite));
    Partial_gw_wrt_theta1 = temp_mat * gw1w2_matrix;
    
    % Get the partial derivative of f with respect to theta_{2}
    temp_mat = diag((0:1:m2)./theta2(ite));
    Partial_fw_wrt_theta2 =  fw1w2_matrix * temp_mat;
    
    % Get the partial derivative of g with respect to theta_{2}
    temp_mat = diag((0:1:n2)./theta2(ite));
    Partial_gw_wrt_theta2 =  gw1w2_matrix * temp_mat;
    
    % Calculate the Partial derivative of T with respect to alpha.
    T1_wrt_alpha = BuildT1(Partial_fw_wrt_alpha,n1,n2,t1,t2,1,1);
    T2_wrt_alpha = BuildT1(Partial_alpha_gw_wrt_alpha,m1,m2,t1,t2,1,1);
    Partial_DTQ_wrt_alpha = D*[T1_wrt_alpha T2_wrt_alpha]*Q;
    
    % Calculate the partial derivative of DTQ with respect to theta_{1}
    T1_wrt_theta1 = BuildT1(Partial_fw_wrt_theta1,n1,n2,t1,t2,1,1);
    T2_wrt_theta1 = BuildT1(Partial_gw_wrt_theta1,m1,m2,t1,t2,1,1);
    Partial_DTQ_wrt_theta1 = D*[T1_wrt_theta1 alpha(ite)*T2_wrt_theta1]*Q;
    
    % Calculate the partial derivative of DTQ with respect to theta_{2}
    T1_wrt_theta2 = BuildT1(Partial_fw_wrt_theta2,n1,n2,t1,t2,1,1);
    T2_wrt_theta2 = BuildT1(Partial_gw_wrt_theta2,m1,m2,t1,t2,1,1);
    Partial_DTQ_wrt_theta2 = D*[T1_wrt_theta2 alpha(ite)*T2_wrt_theta2]*Q;
    
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ck = DTQ*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha, \theta_{1}
    % and \theta_{2}
    Partial_ck_wrt_alpha        = Partial_DTQ_wrt_alpha*e;
    Partial_ck_wrt_theta1       = Partial_DTQ_wrt_theta1*e;
    Partial_ck_wrt_theta2       = Partial_DTQ_wrt_theta2*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fx      = zk(1:num_coeff_f);
    z_gx      = zk(num_coeff_f + 1 :end);
    
    % Get the vectors z_fx and z_gx as matrices, which match the shape of
    % f(x) and g(x).
    z_fx_mat = getAsMatrix(z_fx,m1,m2);
    z_gx_mat = getAsMatrix(z_gx,n1,n2);
    
    % Get matrices z_fw_mat and z_gw_mat, by multiplying rows by
    % theta_{1}^{i} and columns by theta_{2}^{j}
    z_fw_mat = diag(theta1(ite).^(0:1:m1)) * z_fx_mat * diag(theta2(ite).^(0:1:m2));
    z_gw_mat = diag(theta1(ite).^(0:1:n1)) * z_gx_mat * diag(theta2(ite).^(0:1:n2));
    
    % Calculate the derivatives of z_fw and z_gw with repect to \alpha.
    Partial_zfw_wrt_alpha    = zeros(m1+1,m2+1);
    Partial_zgw_wrt_alpha    = z_gw_mat;
    
    % Calculate the derivative of z_fw with respect to \theta_{1}.
    temp_mat = diag((0:1:m1)./theta1(ite));
    Partial_zfw_wrt_theta1 = temp_mat * z_fw_mat;
    
    % Calculate the derivative of z_fw with respect to \theta_{2}
    temp_mat = diag((0:1:m2)./theta2(ite));
    Partial_zfw_wrt_theta2 = z_fw_mat * temp_mat;
    
    % Calculate the derivative of z_gw with respect ot theta1
    temp_mat = diag((0:1:n1)./theta2(ite));
    Partial_zgw_wrt_theta1 = temp_mat * z_gw_mat;
    
    % Calculate the deriviate of z_gw with respect to theta2
    temp_mat = diag((0:1:n2)./theta2(ite));
    Partial_zgw_wrt_theta2 = z_gw_mat * temp_mat;
    
    
    % Build the coefficient Matrix N = [T(z1) T(z2)], of structured perturbations, with
    % same structure as DTQ.
    N1 = BuildT1(z_fw_mat,n1,n2,t1,t2,1,1);
    N2 = BuildT1(z_gw_mat,m1,m2,t1,t2,1,1);
    DNQ = [N1 alpha(ite).*N2];
    
    % Build the coefficient matrix N with respect to alpha
    N1_wrt_alpha = BuildT1(Partial_zfw_wrt_alpha, n1,n2,t1,t2,1,1);
    N2_wrt_alpha = BuildT1(Partial_zgw_wrt_alpha, m1,m2,t1,t2,1,1);
    Partial_DNQ_wrt_alpha = D*[N1_wrt_alpha N2_wrt_alpha]*Q;
    
    % Calculate the derivatives of DNQ with respect to theta
    N1_wrt_theta1 = BuildT1(Partial_zfw_wrt_theta1,n1,n2,t1,t2,1,1);
    N2_wrt_theta1 = BuildT1(Partial_zgw_wrt_theta1,m1,m2,t1,t2,1,1);
    Partial_DNQ_wrt_theta1 = D*[N1_wrt_theta1 alpha(ite).*N2_wrt_theta1]*Q;
    
    % Calculate the derivatives of DNQ with respect to theta
    N1_wrt_theta2 = BuildT1(Partial_zfw_wrt_theta2,n1,n2,t1,t2,1,1);
    N2_wrt_theta2 = BuildT1(Partial_zgw_wrt_theta2,m1,m2,t1,t2,1,1);
    Partial_DNQ_wrt_theta2 = D*[N1_wrt_theta2 alpha(ite).*N2_wrt_theta2]*Q;
    
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{k} the column of S_{k} moved to the RHS
    h = DNQ*e;
    
    % Calculate the derivative of h with respect to alpha
    h_alpha = Partial_DNQ_wrt_alpha*e;
    
    % Calculate the derivative of h with respect to theta1
    h_theta1 = Partial_DNQ_wrt_theta1*e;
    
    % Calculate the derivative of h with respect to theta2
    h_theta2 = Partial_DNQ_wrt_theta2*e;
    
    % Build the matrix (T+N)  
    TN1 = BuildT1(fw1w2_matrix + z_fw_mat,n1,n2,t1,t2,1,1);
    TN2 = BuildT1(gw1w2_matrix + z_gw_mat,m1,m2,t1,t2,1,1);
    DTNQ = D*[TN1 alpha(ite).*TN2]*Q;
    
    % Calculate the paritial derivative of (T+N) with respect to
    % alpha
    TN1_wrt_alpha = BuildT1(Partial_fw_wrt_alpha + Partial_zfw_wrt_alpha, n1,n2,t1,t2,1,1);
    TN2_wrt_alpha = BuildT1(Partial_alpha_gw_wrt_alpha + Partial_zgw_wrt_alpha, m1,m2,t1,t2,1,1);
    DTNQ_alpha = D*[TN1_wrt_alpha TN2_wrt_alpha]*Q;
    
    
    % Calculate the paritial derivative of (T+N) with respect to theta1
    TN1_wrt_theta1 = BuildT1(Partial_fw_wrt_theta1 + Partial_zfw_wrt_theta1,n1,n2,t1,t2,1,1);
    TN2_wrt_theta1 = BuildT1(Partial_gw_wrt_theta1 + Partial_zgw_wrt_theta1,m1,m2,t1,t2,1,1);
    DTNQ_theta1 = D*[TN1_wrt_theta1 alpha(ite).*TN2_wrt_theta1]*Q;
    
    % Calculate the paritial derivative of (T+N) with respect to theta2    
    TN1_wrt_theta2 = BuildT1(Partial_fw_wrt_theta2 + Partial_zfw_wrt_theta2,n1,n2,t1,t2,1,1);
    TN2_wrt_theta2 = BuildT1(Partial_gw_wrt_theta2 + Partial_zgw_wrt_theta2,m1,m2,t1,t2,1,1);
    DTNQ_theta2 = D*[TN1_wrt_theta2 alpha(ite).*TN2_wrt_theta2]*Q;
    
    
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    Y = BuildY(m1,m2,n1,n2,t1,t2,opt_col,xk,alpha(ite),theta1(ite),theta2(ite));
    
    DYG = D*Y*G; 
    
    test1 = DYG * [fxy_vec;gxy_vec];
    first_part = xk(1:(opt_col-1));
    second_part = xk(opt_col:end);
    x = [first_part ; 0 ; second_part];
    test2 = DTQ * x;
    
    
    % Calculate the matrix DP where P is the matrix such that c = P[f;g]
    P = BuildP(m1,m2,n1,n2,theta1(ite),theta2(ite),alpha(ite),t1,t2,opt_col);
    P = D*P*G;
    
    % Get residual as a vector
    rk = (ck+h) - DTNQ*M*xk ;

    % Create the matrix C. This is made up of five submatrices, HZ, Hx,
    % H_alpha and H_theta1 and H_theta2.
    
    Hz          = DYG-P;
    
    Hx          = DTNQ*M;
    
    H_alpha     = DTNQ_alpha*M*xk - (Partial_ck_wrt_alpha + h_alpha);
    
    H_theta1    = DTNQ_theta1*M*xk - (Partial_ck_wrt_theta1 + h_theta1);
    
    H_theta2    = DTNQ_theta2*M*xk - (Partial_ck_wrt_theta2 + h_theta2);
    
    C = [Hz,Hx,H_alpha,H_theta1, H_theta2];  % the matrix C
    
    % Calculate the new right hand vector
    ek = ck + h;
    
    % Update gnew - used in lse problem
    res_vec = rk;
        
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(rk) / norm(ek);

    % Update fnew - used in LSE Problem.
    p = -(yy-start_point);

    
end

% Plot Graphs
figure('name','Residuals in SNTLN')
hold on
title('Residuals in SNTLN')
xlabel('Iterations')
ylabel('log_{10} Residuals')
plot((1:1:ite),log10(condition),'-s')
hold off

figure('name','Theta Variation over Newton Raphson Iterations')
hold on
title('Variation of \theta over Newton Raphson Iteration')
plot((1:1:ite),log10(theta1),'-s')
plot((1:1:ite),log10(theta2),'-s')
hold off

figure('name','Alpha variation over Newton Raphson Iterations')
hold on
title('Alpha variation over Newton Raphson Iterations')
plot((1:1:ite),log10(alpha),'-s')
hold off


%%
if ite == max_iterations
    fprintf('SNTLN Failed to converge, default to input values\n')
    fxy_output = fxy_matrix_n;
    gxy_output = gxy_matrix_n;
    alpha_output = initial_alpha;
    theta1_output = initial_theta_1;
    theta2_output = initial_theta_2;
    X_output = x_ls;
    return;
end

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.

% get the vector zk

zPert_f_vec = zk(1:num_coeff_f);
zPert_f_mat = getAsMatrix(zPert_f_vec,m1,m2);

zPert_g_vec = zk(num_coeff_f+1:end);
zPert_g_mat = getAsMatrix(zPert_g_vec,n1,n2);

fxy_output = fxy_matrix_n + zPert_f_mat;

gxy_output = gxy_matrix_n + zPert_g_mat;

X_output  = xk;

alpha_output = alpha(ite);

theta1_output = theta1(ite);

theta2_output = theta2(ite);


% Print the number of iterations
fprintf('--------------------------------------------------------------------------- \n')
fprintf('Iterations over Sylvester Matrix : %i \n', ite);
fprintf('--------------------------------------------------------------------------- \n')
end

















