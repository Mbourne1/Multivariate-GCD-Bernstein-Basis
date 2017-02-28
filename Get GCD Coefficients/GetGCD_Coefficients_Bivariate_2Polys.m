function [dxy] = GetGCD_Coefficients_Bivariate_2Polys(uxy, vxy, fxy, gxy, t1, t2)
% GetGCD_Coefficients(uxy_matrix,vxy_matrix,fxy_matrix,gxy_matrix,t1,t2,...
%   lambda,mu,alpha, th1,th2)
%
% Given the two quotient polynomials u(x,y) and v(x,y) calculate the GCD d(x,y)
%
% % Inputs 
%
% [uxy, vxy] : (Matrix) Coefficients of polynomial u(x,y) and v(x,y)
%
% [fxy, gxy] : (Matrix) Coefficients of polynomial f(x,y) and g(x,y)
%
% t1 : (Int) Degree of d(x,y) with respect to x
%
% t2 : (Int) Degree of d(x,y) with respect to y
%
% % Outputs
%
% dxy : (Matrix) Coefficients of polynomial d(x,y)


% Get the degrees of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get the degrees of polynomial g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% % Assemble the Subresultant S_{t1,t2}

% Build matrix H
H = BuildH_Bivariate(m1, m2, n1, n2);

% Build the matrix C_{}(u(x,y))
C1_u = BuildT1_Bivariate(uxy, t1, t2);

% Build Matrix C_{}(v(x,y)
C2_v = BuildT1_Bivariate(vxy, t1, t2);

C = [ C1_u;
    C2_v ];

% Buid matrix G
G = BuildG_Bivariate(t1, t2);

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


%
% Arrange dw into a matrix form based on its dimensions.
dxy = GetAsMatrix(x_ls, t1, t2);




end

