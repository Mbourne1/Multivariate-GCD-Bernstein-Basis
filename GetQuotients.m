function [uxy_matrix_calc, vxy_matrix_calc,...
    lambda,mu,...
    opt_alpha, opt_theta_1, opt_theta_2] = GetQuotients(fxy_matrix, gxy_matrix,...
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
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get the degrees of polynomial g(x,y)
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;


%%
switch BOOL_PREPROC
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
    otherwise
        error('err')
end

% Build the (t1,t2)-th subresultant
St1t2 = BuildSubresultant(fxy_matrix_n,gxy_matrix_n,t1,t2,opt_alpha, opt_theta_1, opt_theta_2);


%% Find Optimal column for removal from St
% given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized

[~,cols] = size(St1t2);

% QR Decomposition of the Sylvester Matrix S_{k}
[Qk,Rk] = qr(St1t2);
n = n1+n2;
for j=1:1:cols
    
    ck = St1t2(:,j);
    [Q,~] = qrdelete(Qk,Rk,j);
    cd = Q'*ck;
    d = cd(n+1:end,:);
    residuals_QR(j) = norm(d);
    
end

%Obtain the column for which the residual is minimal.
[~,opt_col] = min(log10(residuals_QR));
fprintf('Optimal column for removal is given by %i \n',opt_col)


%% Get the coefficients for u(x,y) and v(x,y)

Atj = St1t2;
cki = St1t2(:,opt_col);
Atj(:,opt_col) = [];

% [~,n_col] = size(Atj);
% [Q,R] = qr(Atj);
% R1 = R(1:n_col,:);
% cd = Q'*cki;
% c = cd(1:n_col,:);
% x_ls = R1\c;

x_ls = pinv(Atj) * cki;

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ];

num_coeff_v = (n1-t1+1) * (n2-t2+1);
num_coeff_u = (m1-t1+1) * (m2-t2+1);


% get coefficients of u and v
vw_calc = vecx(1:num_coeff_v);
uw_calc = -vecx(num_coeff_v+1:end);


%% Get the value of fv-gu
norm(St1t2 * [vw_calc ;-uw_calc])



%% Obtain u(x,y) in its matrix form
% Arrange uw into a matrix form based on its dimensions.
uw_calc_mat = GetAsMatrix(uw_calc,m1-t1,m2-t2);

%% Obtain v(x,y) in its matrix form
% Arrange vw into a matrix form based on their dimensions.
vw_calc_mat = GetAsMatrix(vw_calc,n1-t1,n2-t2);

%%
% Remove the thetas from the matrix of coefficients of v(w,w) to obtain
% coefficients of v(x,y)

pre_theta = diag(1./(opt_theta_1.^(0:1:n1-t1)));
post_theta = diag(1./(opt_theta_2.^(0:1:n2-t2)));

vxy_matrix_calc = pre_theta * vw_calc_mat * post_theta;

%%
% Remove the thetas from the matrix of coefficients of u(w,w) to obtain
% coefficients of u(x,y)

% for each row, divide by theta2^i1

pre_theta = diag(1./(opt_theta_1.^(0:1:m1-t1)));
post_theta = diag(1./(opt_theta_2.^(0:1:m2-t2)));

uxy_matrix_calc = pre_theta * uw_calc_mat * post_theta;

%% If we excluded Q from the coefficient matrix, then remove the binomial 
% coefficients from v(x,y) and u(x,y)
switch BOOL_Q
    case 'y'
    case 'n'
                %%
        % Remove binomial coefficients from v(w,w)_bi
        [r,c] = size(vxy_matrix_calc);
        n1_t1 = r - 1;
        n2_t2 = c - 1;
        
        bi_n1_t1 = GetBinomials(n1_t1)
        mat1 = diag(1./bi_n1_t1);
        
        bi_n2_t2 = GetBinomials(n2_t2)
        mat2 = diag(1./bi_n2_t2);
        
        vxy_matrix_calc = mat1 * vxy_matrix_calc * mat2;
        
        %%
        % Remove the binomial coefficients from u(w,w)_bi
        [r,c] = size(uxy_matrix_calc);
        m1_t1 = r - 1;
        m2_t2 = c - 1;
        
        bi_m1_t1 = GetBinomials(m1_t1)
        mat1 = diag(1./bi_m1_t1);
        
        bi_m2_t2 = GetBinomials(m2_t2)
        mat2 = diag(1./bi_m2_t2);
        
        uxy_matrix_calc = mat1 * uxy_matrix_calc * mat2;

    otherwise
        error('err')
end

end