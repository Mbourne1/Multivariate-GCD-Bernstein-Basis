function [dxy] = GetGCD_Coefficients_Bivariate_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, t1, t2)
% GetGCD_Coefficients(uxy_matrix,vxy_matrix,fxy_matrix,gxy_matrix,t1,t2,...
%   lambda,mu,alpha, th1,th2)
%
% Given the two quotient polynomials u(x,y) and v(x,y) calculate the GCD d(x,y)
%
% % Inputs 
%
% [uxy, vxy, wxy] : Coefficients of polynomial u(x,y), v(x,y) and w(x,y)
%
% [fxy, gxy, hxy] : Coefficients of polynomial f(x,y), g(x,y) and h(x,y)
%
% k1 : Degree of d(x,y) with respect to x
%
% k2 : Degree of d(x,y) with respect to y
%
% % Outputs
%
% dxy : Coefficients of polynomial d(x,y)


% Get the degrees of polynomial f(x,y), g(x,y) and h(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);


% % Assemble the Subresultant S_{t1,t2}

% Build matrix H
H = BuildH_Bivariate_3Polys(m1, m2, n1, n2, o1, o2);

% Build C_{1}(u)
C1_u = BuildT1_Bivariate(uxy, t1, t2);

% Build Matrix C_{2}(v)
C2_v = BuildT1_Bivariate(vxy, t1, t2);

% Build matrix C_{}(w)
C3_w = BuildT1_Bivariate(wxy, t1, t2);

C = [... 
        C1_u;
        C2_v;
        C3_w;
    ];

% Buid matrix G
G = BuildG_Bivariate(t1, t2);

% Build the Coefficient Matrix HCG 
HCG = H*C*G;

% Get f(x,y), g(x,y) and h(x,y) as vectors
fxy_vec = GetAsVector(fxy);
gxy_vec = GetAsVector(gxy);
hxy_vec = GetAsVector(hxy);

% Create the right hand side vector
rhs_vec = [...
    fxy_vec;
    gxy_vec;
    hxy_vec];


% Obtain vector x
x_ls = SolveAx_b(HCG, rhs_vec);

%
% Arrange dw into a matrix form based on its dimensions.
dxy = GetAsMatrix(x_ls, t1, t2);




end

