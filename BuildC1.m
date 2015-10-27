function C1 = BuildC1(uxy_matrix,t1,t2,m1,m2,th1,th2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                         Inputs


% uxy_matrix :  Input quotient polynomial uxy in matrix form. Excluding
%               binomial coefficients. Excluding thetas.

% t1 :  Degree of GCD with respect to x

% t2 :  Degree of GCD with respect to y

% m1 :  Degree of polynomial f with respect to x

% m2 :  Degree of polynomial f with respect to y

% th1 : Optimal value of theta x = \theta_1\omega_1

% th2 : Optimal value of theta y = \theta_2\omega_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get size of uxy_matrix
[rows,cols] = size(uxy_matrix);

% get m1-t1, degree of uxy with respect to x
m1_t1 = rows - 1;

% Get m2-t2, degree of uxy with respect to y
m2_t2 = cols - 1;

%% Build the matrix uxy_matrix_binoms

% Build the matrix which premultiplies u_xy
Pre_binoms = zeros(m1_t1+1,1);
for i = 0:1:m1_t1
    Pre_binoms(i+1) = nchoosek(m1_t1,i);
end
Pre_binoms_mtrx = diag(Pre_binoms);

% Build the matrix which post multiplies u_xy
Post_binoms = zeros(m2_t2,1);
for i = 0:1:m2_t2
    Post_binoms(i+1) = nchoosek(m2_t2,i);
end
Post_binoms_mtrx = diag(Post_binoms);

% Build the matrix which premultiplies u(x,y) by thetas
Pre_thetas = zeros(m1_t1,1);
for i = 0:1:m1_t1
    Pre_thetas(i+1) = th1^i;
end
Pre_thetas_mtrx = diag(Pre_thetas);

% Build the matrix which post multiplies u(x,y) by thetas
Post_thetas = zeros(m2_t2,1);
for i = 0:1:m2_t2
    Post_thetas(i+1) = th2^i;
end
Post_thetas_mtrx = diag(Post_thetas);

% Create matrix uxy_binom which includes binomial coefficients
uxy_matrix_bi = Pre_binoms_mtrx * uxy_matrix * Post_binoms_mtrx;

% Include thetas
uw1w2_matrix_bi = Pre_thetas_mtrx * uxy_matrix_bi * Post_thetas_mtrx;

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
            temp_vec = getAsVector(uxy_matrix_bi_padded);
            
            C1 = [C1 temp_vec];
        end
    end
    
end