function [dxy_calc_matrix] = GetGCD_Coefficients(uxy_matrix,vxy_matrix,...
    fxy_matrix,gxy_matrix,...
    t1,t2,...
    lambda,mu,...
    opt_alpha, opt_theta_1,opt_theta_2)
%% Given the two quotient polynomials u(x,y) and v(x,y) calculate the GCD d(x,y)

% Inputs 

% uxy_matrix  :
%
% vxy_matrix  :
%
% fxy_matrix  :
%
% gxy_matrix  :
%
% t1          :
%
% t2          :
%
% lambda      :
%
% mu          :
%
% opt_alpha   :
%
% opt_theta_1 :
%
% opt_theta_2 :

%%
% Get the degrees of polynomial f
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get the degrees of polynomial g
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

%% Assemble the Subresultant S_{t1,t2}

% Build matrix H
H = BuildH(m1,m2,n1,n2);

%% Build Matrix HCG
% Build Matrix C1
C1_u = BuildC1(uxy_matrix,t1,t2,m1,m2,opt_theta_1,opt_theta_2);

% Build Matrix C2
C2_v = BuildC1(vxy_matrix,t1,t2,n1,n2,opt_theta_1,opt_theta_2);

C = [ C1_u;
    C2_v ];

% Buid matrix G
G = BuildG(t1,t2);

% Build the Coefficient Matrix HCG 
HCG = H*C*G;

%% Preprocess f(x,y) and g(x,y)

% Noramlise f(x,y) by geometric mean
fxy_matrix_n = fxy_matrix./lambda;

% Noramlise g(x,y) by geometric mean
gxy_matrix_n = gxy_matrix./mu;

% Build the theta matrices which will convert f(x,y) to f(w,w)
pre_theta_mat = diag(opt_theta_1.^(0:1:m1));
post_theta_mat = diag(opt_theta_2.^(0:1:m2));
fww_matrix = pre_theta_mat * fxy_matrix_n * post_theta_mat;

% Get f(x,y) as a vector
fww_vec_th = GetAsVector(fww_matrix);

%% Include thetas in g(x,y) to obtain g(w,w)

% Build the theta matrices which will convert g(x,y) to g(w,w)
pre_theta_mat   = diag(opt_theta_1.^(0:1:n1));
post_theta_mat  = diag(opt_theta_2.^(0:1:n2));
gww_matrix_th = pre_theta_mat * gxy_matrix_n * post_theta_mat;

% Get g(x,y) as a vector
gww_vec = GetAsVector(gww_matrix_th);


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
dww_calc_matrix = GetAsMatrix(dww_calc,t1,t2);

%%
% Remove thetas from dw
% Divide the row i by theta1^i
mat1 = diag(1./opt_theta_1.^(0:1:t1));
mat2 = diag(1./opt_theta_2.^(0:1:t2));

dxy_calc_matrix = mat1 * dww_calc_matrix * mat2;




end

