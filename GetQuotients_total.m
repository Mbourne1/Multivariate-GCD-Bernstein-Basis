function [uxy_matrix_calc, vxy_matrix_calc] = GetQuotients_total(fxy_matrix, gxy_matrix,...
    m,n,t)
% Given two input polynomials and the degree of the GCD return the quotient
% polynomials u(x,y) and v(x,y)
%
%       Inputs.
%
%   fxy_matrix :
%
%   gxy_matrix :
%
%   m :
%
%   n :
%
%   t :
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialise Global Variables
global bool_Q
global bool_preproc


% padd polynomials fxy and gxy
[m1,m2] = GetDegree(fxy_matrix);
fxy_padd = zeros(m+1,m+1);
fxy_padd(1:m1+1,1:m2+1) = fxy_matrix;

[n1,n2] = size(gxy_matrix);
gxy_padd = zeros(n+1,n+1);
gxy_padd(1:n1,1:n2) = gxy_matrix;

fxy_matrix = fxy_padd;
gxy_matrix = gxy_padd;



%%
switch bool_preproc
    case 'y'
        
        % Preproecessor One - Normalise by geometric mean
        [lambda, mu] = getGeometricMean(fxy_matrix,gxy_matrix,t1,t2);
        
        
        % Normalise f(x,y)
        fxy_matrix_n = fxy_matrix./lambda;
        
        % Normalise g(x,y)
        gxy_matrix_n = gxy_matrix./mu;
        
        % Preprocessor Two and Three - LinProg to obtain optimal values
        % of alpha, theta_1 and theta_2
        
        % Get the maximum and minimum entries of f(x,y) in the
        % Sylvester matrix S(f,g)
        [fxy_max_mtrx,fxy_min_mtrx] = GetMaxMin(fxy_matrix_n,n1,n2,t1,t2);
        
        % Get the maximum and minimum entries of g(x,y) in the
        % Sylvester matrix S(f,g)
        [gxy_max_mtrx,gxy_min_mtrx] = GetMaxMin(gxy_matrix_n,m1,m2,t1,t2);
        
        % Get optimal values of alpha and theta
        [opt_alpha, opt_theta_1, opt_theta_2] = ...
            OptimalAlphaTheta(fxy_max_mtrx,fxy_min_mtrx,gxy_max_mtrx,gxy_min_mtrx);
        

    case 'n'
        opt_theta_1 = 1;
        opt_theta_2 = 1;
        opt_alpha = 1;
        
        
        lambda = 1;
        mu = 1;
        
        fxy_matrix_n = fxy_matrix./lambda;
        gxy_matrix_n = gxy_matrix./mu;
        
end

% Build the (t1,t2)-th subresultant
St1t2 = BuildSubresultant(fxy_matrix_n,gxy_matrix_n,t,t,opt_alpha, opt_theta_1, opt_theta_2);


%%

[opt_col] = GetOptimalColumn(fxy_matrix,gxy_matrix,t1,t2,lambda,mu,opt_alpha,opt_theta_1,opt_theta_2);

Atj = St1t2;
cki = St1t2(:,opt_col);
Atj(:,opt_col) = [];

x_ls = SolveAx_b(Atj,cki);

%% Get the coefficients for u(x,y) and v(x,y)

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ];

num_coeff_v = (n-t+1) * (n-t+1);

% get coefficients of u and v
vw_calc = vecx(1:num_coeff_v);
uw_calc = -vecx(num_coeff_v+1:end);


%% Obtain u(x,y) in its matrix form
% Arrange uw into a matrix form based on its dimensions.
uw_calc_mat = GetAsMatrix(uw_calc,m-t,m-t);

%% Obtain v(x,y) in its matrix form
% Arrange vw into a matrix form based on their dimensions.
vw_calc_mat = GetAsMatrix(vw_calc,n-t,n-t);


% Remove the thetas from the matrix of coefficients of v(w,w) to obtain
% coefficients of v(x,y)
vxy_matrix_calc = GetWithoutThetas(vw_calc_mat,opt_theta_1,opt_theta_2);

% Remove the thetas from the matrix of coefficients of u(w,w) to obtain
% coefficients of u(x,y)
uxy_matrix_calc = GetWithoutThetas(uw_calc_mat,opt_theta_1,opt_theta_2);


% If we excluded Q from the coefficient matrix, then remove the binomial 
% coefficients from v(x,y) and u(x,y)
switch bool_Q
    case 'n'
        
        % Remove binomial coefficients from v(w,w)_bi
        vxy_matrix_calc = GetWithoutBinomials(vxy_matrix_calc);
       
        % Remove the binomial coefficients from u(w,w)_bi
        uxy_matrix_calc = GetWithoutBinomials(uxy_matrix_calc);
       

end

end