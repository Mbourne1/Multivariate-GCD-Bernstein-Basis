function [dxy_calc_matrix] = GetGCD_Coefficients(uxy_matrix,vxy_matrix,...
    fxy_matrix,gxy_matrix,...
    m,n,t,...
    lambda,mu,...
    opt_alpha, opt_theta_1,opt_theta_2)
%% Given the two quotient polynomials u(x,y) and v(x,y) calculate the GCD d(x,y)

% Inputs 

% uxy_matrix  :
% vxy_matrix  :
% fxy_matrix  :
% gxy_matrix  :
% t1          :
% t2          :
% lambda      :
% mu          :
% opt_alpha   :
% opt_theta_1 :
% opt_theta_2 :

%%
% Get the degrees of polynomial f
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

fxy_padd = zeros(m+1,m+1)
fxy_padd(1:m1+1,1:m2+1) = fxy_matrix;
fxy_matrix = fxy_padd;


% Get the degrees of polynomial g
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

gxy_padd = zeros(m+1,m+1)
gxy_padd(1:n1+1,1:n2+1) = gxy_matrix;
gxy_matrix = gxy_padd;


%% Assemble the Subresultant S_{t1,t2}

% Build matrix H
H = BuildH(m,m,n,n);

% % Build Matrix C
% Build Matrix C1
C1_u = BuildC1(uxy_matrix,t,t,m,m,opt_theta_1,opt_theta_2);

% Build Matrix C2
C2_v = BuildC1(vxy_matrix,t,t,n,n,opt_theta_1,opt_theta_2);

C = [ C1_u;
    C2_v ];

% Buid matrix G
G = BuildG(t,t);

% Build the Coefficient Matrix HCG 
HCG = H*C*G;

%% Preprocess f(x,y) and g(x,y)

% Noramlise f(x,y) by geometric mean
fxy_matrix_n = fxy_matrix./lambda;

% Noramlise g(x,y) by geometric mean
gxy_matrix_n = gxy_matrix./mu;

% Build the theta matrices which will convert f(x,y) to f(w,w)
pre_theta_mat = diag(opt_theta_1.^(0:1:m));
post_theta_mat = diag(opt_theta_2.^(0:1:m));
fww_matrix = pre_theta_mat * fxy_matrix_n * post_theta_mat;

% Get f(x,y) as a vector
fww_vec_th = getAsVector(fww_matrix);

%% Include thetas in g(x,y) to obtain g(w,w)

% Build the theta matrices which will convert g(x,y) to g(w,w)
pre_theta_mat   = diag(opt_theta_1.^(0:1:n));
post_theta_mat  = diag(opt_theta_2.^(0:1:n));
gww_matrix_th = pre_theta_mat * gxy_matrix_n * post_theta_mat;

% Get g(x,y) as a vector
gww_vec = getAsVector(gww_matrix_th);


%% Create the right hand side vector
rhs_vec = [...
    fww_vec_th;
    opt_alpha.*gww_vec];


%% Obtain vector x
x = pinv(HCG) * rhs_vec;
dww_calc = x;

residual = pinv(HCG)*rhs_vec - x;

%%
% Arrange dw into a matrix form based on its dimensions.
dww_calc_matrix = getAsMatrix(dww_calc,t,t);

%%
% Remove thetas from dw
% Divide the row i by theta1^i
mat1 = diag(1./opt_theta_1.^(0:1:t));
mat2 = diag(1./opt_theta_2.^(0:1:t));

dxy_calc_matrix = mat1 * dww_calc_matrix * mat2;




end

