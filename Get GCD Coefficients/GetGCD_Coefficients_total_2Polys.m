function [dxy_calc_matrix] = GetGCD_Coefficients_total_2Polys(uxy,vxy,...
    fxy,gxy,...
    m,n,t,...
    GM_fx,GM_gx,...
    alpha, th1,th2)
% Given the two quotient polynomials u(x,y) and v(x,y) calculate the GCD d(x,y)
%
% Inputs 
%
% [uxy, vxy] :
%
% [fxy, gxy] :
%
% t1 :
%
% t2 :
%
% GM_fx :
%
% GM_gx :
%
% opt_alpha :
%
% opt_theta_1 :
%
% opt_theta_2 :

%
% Get the degrees of polynomial f
[m1,m2] = GetDegree(fxy);

fxy_padd = zeros(m+1,m+1);
fxy_padd(1:m1+1,1:m2+1) = fxy;
fxy = fxy_padd;


% Get the degrees of polynomial g
[n1,n2] = GetDegree(gxy);

gxy_padd = zeros(m+1,m+1);
gxy_padd(1:n1+1,1:n2+1) = gxy;
gxy = gxy_padd;

% Assemble the Subresultant S_{t1,t2}

% Build matrix H
H = BuildH(m,m,n,n);

% % Build Matrix C
% Build Matrix C1
uww_matrix = GetWithThetas(uxy,th1,th2);
vww_matrix = GetWithThetas(vxy,th1,th2);

C1_u = BuildT1(uww_matrix,t,t);
C2_v = BuildT1(vww_matrix,t,t);

C = [ ...
    C1_u;
    C2_v;
    ];

% Buid matrix G
G = BuildG(t,t);

% Build the Coefficient Matrix HCG 
HCG = H*C*G;

% %
% Preprocess f(x,y) and g(x,y)

% Noramlise f(x,y) by geometric mean
fxy_matrix_n = fxy./GM_fx;

% Noramlise g(x,y) by geometric mean
gxy_matrix_n = gxy./GM_gx;

% Build the theta matrices which will convert f(x,y) to f(w,w)
fww_matrix = GetWithThetas(fxy_matrix_n);

% Get f(w,w) as a vector
fww_vec_th = GetAsVector(fww_matrix);

% % Include thetas in g(x,y) to obtain g(w,w)
%

% Build the theta matrices which will convert g(x,y) to g(w,w)
gww_matrix_th = GetWithThetas(gxy_matrix_n,th1,th2);

% Get g(x,y) as a vector
gww_vec = GetAsVector(gww_matrix_th);


%% Create the right hand side vector
rhs_vec = [...
    fww_vec_th;
    alpha.*gww_vec];


%% Obtain vector x

x_ls = SolveAx_b(HCG,rhs_vec);

dww_calc = x_ls;

%residual = pinv(HCG)*rhs_vec - x;

%%
% Arrange dw into a matrix form based on its dimensions.
dww_calc_matrix = getAsMatrix(dww_calc,t,t);

%%
% Remove thetas from dw
% Divide the row i by theta1^i
dxy_calc_matrix = GetWithoutThetas(dww_calc_matrix);



end

