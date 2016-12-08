function [dxy] = GetGCD_Coefficients_Relative(uxy, vxy, fxy, gxy, t1, t2)
% GetGCD_Coefficients(uxy_matrix,vxy_matrix,fxy_matrix,gxy_matrix,t1,t2,...
%   lambda,mu,alpha, th1,th2)
%
% Given the two quotient polynomials u(x,y) and v(x,y) calculate the GCD d(x,y)
%
% % Inputs 
%
% uxy : Coefficients of polynomial u(x,y)
%
% vxy : Coefficients of polynomial v(x,y)
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% k1 : Degree of d(x,y) with respect to x
%
% k2 : Degree of d(x,y) with respect to y
%
% % Outputs
%
% dxy : Coefficients of polynomial d(x,y)


% Get the degrees of polynomial f
[m1,m2] = GetDegree(fxy);

% Get the degrees of polynomial g
[n1,n2] = GetDegree(gxy);

% % Assemble the Subresultant S_{t1,t2}

% Build matrix H
H = BuildH(m1,m2,n1,n2);

C1_u = BuildT1(uxy,t1,t2);

% Build Matrix C2
C2_v = BuildT1(vxy,t1,t2);

C = [ C1_u;
    C2_v ];

% Buid matrix G
G = BuildG(t1,t2);

% Build the Coefficient Matrix HCG 
HCG = H*C*G;

% Get f(x,y) as a vector
fxy_vec = GetAsVector(fxy);

% Get g(x,y) as a vector
gxy_vec = GetAsVector(gxy);


% Create the right hand side vector
rhs_vec = [...
    fxy_vec;
    gxy_vec];


% Obtain vector x
x_ls = SolveAx_b(HCG,rhs_vec);
dxy_calc = x_ls;

%
% Arrange dw into a matrix form based on its dimensions.
dxy = GetAsMatrix(dxy_calc,t1,t2);




end

