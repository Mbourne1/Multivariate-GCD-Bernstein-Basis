function [dxy] = GetGCD_Coefficients_Bivariate_2Polys(uxy, vxy, fxy, gxy, t1, t2)
% GetGCD_Coefficients(uxy_matrix,vxy_matrix,fxy_matrix,gxy_matrix,t1,t2,...
%   lambda,mu,alpha, th1,th2)
%
% Given the two quotient polynomials u(x,y) and v(x,y) calculate the GCD d(x,y)
%
% % Inputs 
%
% uxy : (Matrix) Coefficients of polynomial u(x,y), the quotient polynomial
% such that f(x,y) / u(x,y) = d(x,y)
%
% vxy : (Matrix) Coefficients of polynomial v(x,y), the quotient polynomial
% such that g(x,y) / v(x,y) = d(x,y)
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) 
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% t1 : (Int) Degree of d(x,y) with respect to x
%
% t2 : (Int) Degree of d(x,y) with respect to y
%
% % Outputs
%
% dxy : (Matrix) Coefficients of polynomial d(x,y), the greatest common
% divisor of f(x,y) and g(x,y)


% % Assemble the factorisation matrix 
C_uv = BuildFactorisationMatrix(uxy, vxy, t1, t2);

% Get f(x,y) as a vector
fxy_vec = GetAsVector(fxy);

% Get g(x,y) as a vector
gxy_vec = GetAsVector(gxy);


% Create the right hand side vector
rhs_vec = [...
    fxy_vec;
    gxy_vec];


% Obtain vector x
x_ls = SolveAx_b(C_uv,rhs_vec);


%
% Arrange dw into a matrix form based on its dimensions.
dxy = GetAsMatrix(x_ls, t1, t2);




end


