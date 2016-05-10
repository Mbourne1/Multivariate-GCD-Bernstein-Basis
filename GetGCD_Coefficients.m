function [dxy_calc_matrix] = GetGCD_Coefficients(uxy_matrix,vxy_matrix,...
    fxy_matrix_n,gxy_matrix_n,...
    t1,t2)
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
% Get the degrees of polynomial f
[m1,m2] = GetDegree(fxy_matrix_n);

% Get the degrees of polynomial g
[n1,n2] = GetDegree(gxy_matrix_n);

% % Assemble the Subresultant S_{t1,t2}

% Build matrix H
H = BuildH(m1,m2,n1,n2);

C1_u = BuildT1(uxy_matrix,t1,t2);

% Build Matrix C2
C2_v = BuildT1(vxy_matrix,t1,t2);

C = [ C1_u;
    C2_v ];

% Buid matrix G
G = BuildG(t1,t2);

% Build the Coefficient Matrix HCG 
HCG = H*C*G;

% % Preprocess f(x,y) and g(x,y)

% Get f(x,y) as a vector
fxy_vec_th = GetAsVector(fxy_matrix_n);

% % Include thetas in g(x,y) to obtain g(w,w)
% Get g(x,y) as a vector
gxy_vec = GetAsVector(gxy_matrix_n);


% Create the right hand side vector
rhs_vec = [...
    fxy_vec_th;
    gxy_vec];


% Obtain vector x
x_ls = SolveAx_b(HCG,rhs_vec);
dxy_calc = x_ls;

%
% Arrange dw into a matrix form based on its dimensions.
dxy_calc_matrix = GetAsMatrix(dxy_calc,t1,t2);




end

