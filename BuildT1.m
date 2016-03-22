function T1 = BuildT1(fxy_matrix,n1_k1,n2_k2)
% Build the matrix T1 where T1(f) * [v] gives a vector of coefficients of
% the polynomial multiplication f(x,y) * v(x,y).
%
% Inputs.
%
% fww_matrix : Coefficients of polynomial f(w,w)
%
% n1_k1 : Degree of polynomial v(x,y) with respect to x
%
% n2_k2 : Degree of polynomial v(x,y) with repsect to y


% get size of f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% % Build the matrix f(w,w) including binomials.
fxy_matrix_bi = GetWithBinomials(fxy_matrix);

% Initalise a zero matrix 
zero_matrix = zeros(m1 + n1_k1 + 1, m2 + n2_k2 + 1);

% Get the number of rows in T_{k_{1},k_{2}}(f)
nColsT1 = (n1_k1 + 1) * (n2_k2 + 1);
nRowsT1 = (m1+n1_k1 +1 ) * (m2+n2_k2 +1);

T1 = zeros(nRowsT1,nColsT1);

% for every diagonal of the matrix vxy. 
nDiags_vxy = (n2_k2+1)+(n1_k1+1) -1;  

% Initialise a count
count = 1;

for tot0 = 0:1:nDiags_vxy - 1;
    for i = tot0:-1:0
        j = tot0-i;
        % if j1 is within the number of columns of vxy_mtrx.
        % if j2 is within the number of rows of vxy_mtrx.
        if j <= n2_k2  && i <= n1_k1 
            
            % Padd the coefficients
            fxy_matrix_bi_padded_new = zero_matrix;
            fxy_matrix_bi_padded_new((i+1):(m1)+(i+1),(j+1):(m2)+(j+1)) = fxy_matrix_bi;
            
            temp_vec = GetAsVector(fxy_matrix_bi_padded_new);
            
            % Assign column
            T1(:,count) = temp_vec;
            
            % Increment Counter
            count = count + 1;
        end
        %
    end
end



end