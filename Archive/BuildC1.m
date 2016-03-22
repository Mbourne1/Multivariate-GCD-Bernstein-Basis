function C1 = BuildC1(uxy_matrix,t1,t2,m1,m2,th1,th2)
% Build the matrix C1(u), where C1(u) * d = f.
%
% %                         Inputs
%
%
% uxy_matrix :  Input quotient polynomial uxy in matrix form. Excluding
%               binomial coefficients. Excluding thetas.
%
% t1 :  Degree of GCD with respect to x
%
% t2 :  Degree of GCD with respect to y
%
% m1 :  Degree of polynomial f with respect to x
%
% m2 :  Degree of polynomial f with respect to y
%
% th1 : Optimal value of theta x = \theta_1\omega_1
%
% th2 : Optimal value of theta y = \theta_2\omega_2
%
% %


% Get size of uxy_matrix
[nRows_uxy,nCols_uxy] = size(uxy_matrix);

% get m1-t1, degree of uxy with respect to x
m1_t1 = nRows_uxy - 1;

% Get m2-t2, degree of uxy with respect to y
m2_t2 = nCols_uxy - 1;

%% Build the matrix uxy_matrix_binoms

%Include Binomial Coefficients
uxy_matrix_bi = GetWithBinomials(uxy_matrix);

% Include thetas
uw1w2_matrix_bi = GetWithThetas(uxy_matrix_bi,th1,th2);

%%

% initialise a zero matrix to pad uxy
zero_matrix = zeros(m1+1,m2+1);

%
C1 = [];

% Get number of diagonals in d(x,y)
num_diags = (t1 + 1) + (t2 + 1) - 1;

% for every diagonal of the matrix dxy_mtrx.
for tot = 0:1:num_diags
    
    for i = tot:-1:0
        j = tot-i;

        % if j1 is within the number of columns of dxy_mtrx
        % if j2 is within the number of rows of dxy_mtrx
        if i <= t1 &&  j <= t2
            
            
            uxy_matrix_bi_padded = zero_matrix;
            
            uxy_matrix_bi_padded((i+1):(m1_t1)+(i+1), (j+1):(m2_t2)+(j+1)) = ...
                uw1w2_matrix_bi;
            
            % Get u(x,y)* basis element, as a vector.
            temp_vec = GetAsVector(uxy_matrix_bi_padded);
            
            C1 = [C1 temp_vec];
        end
    end
    
end