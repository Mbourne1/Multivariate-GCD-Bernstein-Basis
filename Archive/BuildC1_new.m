function T1 = BuildC1_new(uww_matrix_bi,m1,m2,n1,n2,t1,t2)


% Get the degree of u(w,w)
[r,c] = size(uww_matrix_bi);
m1_t1 = r-1;
m2_t2 = c-1;

% Initalise a zero matrix 
zero_matrix = zeros(m1_t1 + n1 + 1, m2_t2 + n2 + 1);


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