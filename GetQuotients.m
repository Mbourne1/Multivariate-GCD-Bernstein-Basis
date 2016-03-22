function [uxy_matrix_calc, vxy_matrix_calc,...
    lambda,mu,...
    opt_alpha, th1, th2] = GetQuotients(fxy_matrix, gxy_matrix,...
    t1,t2)
% Given two input polynomials and the degree of the GCD return the quotient
% polynomials u(x,y) and v(x,y)
%
%   Inputs.
%
%   fxy_matrix :
%
%   gxy_matrix :
%
%   t1 :
%
%   t2 :
%

% Initialise Global Variables
global BOOL_Q

global BOOL_PREPROC

%%
% Get the degrees of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get the degrees of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);



%
switch BOOL_PREPROC
    case 'y'
        
        % Preproecessor One - Normalise by geometric mean
        [lambda, mu] = GetGeometricMean(fxy_matrix,gxy_matrix,t1,t2);
        
        
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
        [opt_alpha, th1, th2] = ...
            OptimalAlphaTheta(fxy_max_mtrx,fxy_min_mtrx,gxy_max_mtrx,gxy_min_mtrx);
        

    case 'n'
        th1 = 1;
        th2 = 1;
        opt_alpha = 1;
        
        
        lambda = 1;
        mu = 1;
        
        fxy_matrix_n = fxy_matrix./lambda;
        gxy_matrix_n = gxy_matrix./mu;
    otherwise
        error('err')
end

% Build the (t1,t2)-th subresultant
St1t2 = BuildSubresultant(fxy_matrix_n,gxy_matrix_n,t1,t2,opt_alpha, th1, th2);


%% Find Optimal column for removal from St
% given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized
opt_col_index = GetOptimalColumn(fxy_matrix,gxy_matrix,t1,t2,lambda,mu,opt_alpha, th1, th2);

% Get the matrix A_{t_{1},t_{2}} 
At = St1t2;
At(:,opt_col_index) = [];

% Get the vector c_{t_{1},t_{2}} removed from S_{t_{1},t_{2}}
ct = St1t2(:,opt_col_index);


%% Get the coefficients for u(x,y) and v(x,y)
x_ls = SolveAx_b(At,ct);


% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col_index)-1);
    -1;
    x_ls(opt_col_index:end);
    ];

num_coeff_v = (n1-t1+1) * (n2-t2+1);
num_coeff_u = (m1-t1+1) * (m2-t2+1);


% get coefficients of u(w,w) and v(w,w)
vww_calc = vecx(1:num_coeff_v);
uww_calc = -vecx(num_coeff_v+1:end);


%% Get the value of fv-gu
norm(St1t2 * [vww_calc ;-uww_calc])

%% Obtain u(x,y) in its matrix form
% Arrange uw into a matrix form based on its dimensions.
uw_calc_mat = GetAsMatrix(uww_calc,m1-t1,m2-t2);

%% Obtain v(x,y) in its matrix form
% Arrange vw into a matrix form based on their dimensions.
vw_calc_mat = GetAsMatrix(vww_calc,n1-t1,n2-t2);

% Remove the thetas from the matrix of coefficients of v(w,w) to obtain
% coefficients of v(x,y)
vxy_matrix_calc = GetWithoutThetas(vw_calc_mat,th1,th2);

% Remove the thetas from the matrix of coefficients of u(w,w) to obtain
% coefficients of u(x,y)
uxy_matrix_calc = GetWithoutThetas(uw_calc_mat,th1,th2);

% If we excluded Q from the coefficient matrix, then remove the binomial 
% coefficients from v(x,y) and u(x,y)

switch BOOL_Q
    case 'y'
        % Do nothing
    case 'n'
                %%
        % Remove binomial coefficients from v(w,w)_bi
        [n1_t1,n2_t2] = GetDegree(vxy_matrix_calc);
        
        vxy_matrix_calc = GetWithoutBinomails(vxy_matrix_calc);
        
        uxy_matrix_calc = GetWithoutBinomials(uxy_matrix_calc);

    otherwise
        error('err')
end

end