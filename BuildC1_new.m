function T1 = BuildC1_new(uww_matrix,m1,m2,n1,n2,t1,t2,~,~)
% % Build the matrix u_xy_binoms by pre and post multiplying 'uxy_matrix'
% by two binomial matrices

% Build the matrix which pre multiplies u(x,y)
Bi_m1_t1 = GetBinomials(m1-t1);
Pre_binoms_mtrx = diag(Bi_m1_t1);

% Build the matrix which post multiplies u(x,y)
Bi_m2_t2 = GetBinomials(m2-t2);
Post_binoms_mtrx = diag(Bi_m2_t2);

% Create matrix fxy_bi which includes the binomial coefficients
uww_matrix_bi = Pre_binoms_mtrx *  uww_matrix * Post_binoms_mtrx;


% Initalise a zero matrix 
zero_matrix = zeros(m1 + n1 - t1 + 1, m2 + n2 - t2 + 1);
uww_matrix_bi_padded = zero_matrix;

uww_matrix_bi_padded(1:m1-t1+1,1:m2-t2+1) = uww_matrix_bi;

T1 = [];

% Get the number of diagonals in g(x,y)
nDiags_gxy = (n1+1) + (n2+1) -1;


% Build the matrix C1(u) by multiplying the coefficients by each of the
% basis elements of g(x,y), in order.

for tot = 0:1: nDiags_gxy - 1
    for i = tot:-1:0
        j = tot-i;
        % if j1 is within the number of columns of vxy_mtrx.
        % if j2 is within the number of rows of vxy_mtrx.
        if j <= n2  && i <= n1 

            uww_matrix_bi_padded_new = zero_matrix;
            
            uww_matrix_bi_padded_new((i+1):(m1-t1)+(i+1),(j+1):(m2-t2)+(j+1)) = uww_matrix_bi;
            
            temp_vec = GetAsVector(uww_matrix_bi_padded_new);
            
            T1 = [T1 temp_vec];
        end
        %
    end
end



end