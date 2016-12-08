function [uxy, vxy] = GetCofactors_Total(fxy, gxy, m, n, t)
% Given two input polynomials and the degree of the GCD return the quotient
% polynomials u(x,y) and v(x,y)
%
% % Inputs.
%
% fxy : Coefficients of the polynomial f(x,y) in the Bernstein basis
%
% gxy : Coefficients of the polynomial g(x,y) in the Bernstein basis
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% t : Total degree of d(x,y)


% Initialise Global Variables
global SETTINGS

% padd polynomials fxy and gxy
[m1,m2] = GetDegree(fxy);
fxy_padd = zeros(m+1,m+1);
fxy_padd(1:m1+1,1:m2+1) = fxy;

[n1,n2] = size(gxy);
gxy_padd = zeros(n+1,n+1);
gxy_padd(1:n1,1:n2) = gxy;

fxy = fxy_padd;
gxy = gxy_padd;


% Build the (t1,t2)-th subresultant
St1t2 = BuildSubresultant(fxy_matrix_n,gxy_matrix_n,t,t,opt_alpha, opt_theta_1, opt_theta_2);


% %

[opt_col] = GetOptimalColumn(fxy,gxy,t1,t2,lambda,mu,opt_alpha,opt_theta_1,opt_theta_2);

Atj = St1t2;
cki = St1t2(:,opt_col);
Atj(:,opt_col) = [];

x_ls = SolveAx_b(Atj,cki);

% % Get the coefficients for u(x,y) and v(x,y)

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ];

nCoeffs_vxy = (n-t+1) * (n-t+1);

% get coefficients of u and v
vw_calc = vecx(1:nCoeffs_vxy);
uw_calc = -vecx(nCoeffs_vxy+1:end);


% % Obtain u(x,y) in its matrix form
% Arrange uw into a matrix form based on its dimensions.
uw_calc_mat = GetAsMatrix(uw_calc,m-t,m-t);

% % Obtain v(x,y) in its matrix form
% Arrange vw into a matrix form based on their dimensions.
vw_calc_mat = GetAsMatrix(vw_calc,n-t,n-t);


% If we excluded Q from the coefficient matrix, then remove the binomial 
% coefficients from v(x,y) and u(x,y)
switch SETTINGS.SYLVESTER_BUILD_METHOD
    case 'T'
        
    case 'DT'
        
        % Remove binomial coefficients from v(w,w)_bi
        vxy = GetWithoutBinomials(vxy);
        
        % Remove the binomial coefficients from u(w,w)_bi
        uxy = GetWithoutBinomials(uxy);
        
    case 'DTQ'
        
    case 'TQ'
        
end

end