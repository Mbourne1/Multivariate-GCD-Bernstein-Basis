
function T1 = BuildT1(fxy_matrix,n1,n2,k1,k2,th1,th2)

% get size of f(x,y)
[rows,cols] = size(fxy_matrix);

% Get degree of f with respect to x
m1 = rows - 1;

% Get degree of f with respect to y
m2 = cols - 1;


% % Build the matrix f_xy_binoms by pre and post multiplying 'fxy_matrix'
% by two binomial matrices

% Build the matrix which pre multiplies f_xy
Pre_binoms = GetBinomials(m1);
Pre_binoms_mtrx = diag(Pre_binoms);

% Build the matrix which post multiplies f_xy
Post_binoms = GetBinomials(m2);
Post_binoms_mtrx = diag(Post_binoms);

% Build the matrix which pre multiplies f_xy
Pre_thetas_mtrx = diag(th1.^(0:1:m1));

% Build the matrix which post multiplies f_xy
Post_thetas_mtrx = diag(th2.^(0:1:m2));

% Create matrix fxy_bi which includes the binomial coefficients
fxy_matrix_bi = Pre_binoms_mtrx *  fxy_matrix * Post_binoms_mtrx;

% Include thetas and binomial coefficients
fw1w2_matrix_bi = Pre_thetas_mtrx * fxy_matrix_bi * Post_thetas_mtrx;

% Initalise a zero matrix 
zero_matrix = zeros(m1 + n1 - k1 + 1, m2 +n2 - k2 + 1);
fxy_matrix_bi_padded = zero_matrix;

T1 = [];

% for every diagonal of the matrix vxy. 
num_diags_vxy = (n2-k2+1)+(n1-k1+1) -1;  

for tot0 = 0:1:num_diags_vxy - 1;
    
    for i = tot0:-1:0
        j = tot0-i;
        % if j1 is within the number of columns of vxy_mtrx.
        % if j2 is within the number of rows of vxy_mtrx.
        if j <= n2-k2  && i <= n1-k1 
            
            fxy_matrix_bi_padded_new = zero_matrix;
            fxy_matrix_bi_padded_new((i+1):(m1)+(i+1),(j+1):(m2)+(j+1)) = fw1w2_matrix_bi;
            
            temp_vec = GetAsVector(fxy_matrix_bi_padded_new);
            
            T1 = [T1 temp_vec];
        end
        %
    end
end

T1 = T1;

end