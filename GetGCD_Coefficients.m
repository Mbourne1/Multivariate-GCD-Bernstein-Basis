function [dxy_calc_matrix] = GetGCD_Coefficients(uxy_matrix,vxy_matrix,...
    fxy_matrix,gxy_matrix,...
    t1,t2,...
    lambda,mu,...
    alpha, th1,th2)
% GetGCD_Coefficients(uxy_matrix,vxy_matrix,fxy_matrix,gxy_matrix,t1,t2,...
%   lambda,mu,alpha, th1,th2)
%
% Given the two quotient polynomials u(x,y) and v(x,y) calculate the GCD d(x,y)
%
% Inputs 
%
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

%
% Get the degrees of polynomial f
[m1,m2] = GetDegree(fxy_matrix);

% Get the degrees of polynomial g
[n1,n2] = GetDegree(gxy_matrix);

% % Assemble the Subresultant S_{t1,t2}

% Build matrix H
H = BuildH(m1,m2,n1,n2);

% % Build Matrix HCG
% Build Matrix C1
uww_matrix = GetWithThetas(uxy_matrix,th1,th2);
vww_matrix = GetWithThetas(vxy_matrix,th1,th2);

C1_u = BuildT1(uww_matrix,t1,t2);

% Build Matrix C2
C2_v = BuildT1(vww_matrix,t1,t2);

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
fww_matrix = GetWithThetas(fxy_matrix_n,th1,th2);

% Get f(x,y) as a vector
fww_vec_th = GetAsVector(fww_matrix);

%% Include thetas in g(x,y) to obtain g(w,w)

% Build the theta matrices which will convert g(x,y) to g(w,w)
gww_matrix = GetWithThetas(gxy_matrix_n,th1,th2);

% Get g(x,y) as a vector
gww_vec = GetAsVector(gww_matrix);


% Create the right hand side vector
rhs_vec = [...
    fww_vec_th;
    alpha.*gww_vec];


% Obtain vector x
x_ls = SolveAx_b(HCG,rhs_vec);
dww_calc = x_ls;

%
% Arrange dw into a matrix form based on its dimensions.
dww_calc_matrix = GetAsMatrix(dww_calc,t1,t2);

%
% Remove thetas from dw
% Divide the row i by theta1^i

dxy_calc_matrix = GetWithoutThetas(dww_calc_matrix,th1,th2);




end

