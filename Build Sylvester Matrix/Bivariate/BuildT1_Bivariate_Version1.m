function T1 = BuildT1_Bivariate_Version1(fxy,n1_k1,n2_k2)
% % Build the matrix T_{n1 - k1, n2 - k2}(f(x,y)), which is the
% (n_{1}-k_{1}, n_{2} - k_{2}) th convolution matrix of the polynomial
% f(x, y)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of f(x,y)
%
% n1_k1 : (Int) Degree of v(x,y) with respect to x
%
% n2_k2 : (Int) Degree of v(x,y) with respect to y


% Get the degree of f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);

% Get coefficients of f(x,y) in scaled Bernstein form
fxy_matrix_bi = GetWithBinomials_Bivariate(fxy);

% Initalise a zero matrix
zero_matrix = zeros(m1 + n1_k1 + 1, m2 + n2_k2 + 1);

% Get the number of rows and columns in T_{k_{1},k_{2}}(f(x,y))
nRowsT1 = (m1 + n1_k1 +1 ) * (m2+n2_k2 +1);
nColumnsT1 = (n1_k1 + 1) * (n2_k2 + 1);

% Initialise a zero matrix T_{n_{1}-k_{1},n_{2}-k_{2}}
T1 = zeros(nRowsT1,nColumnsT1);

% For every diagonal of the matrix vxy.
nDiagonals_vxy = (n2_k2 + 1)+(n1_k1 + 1) -1;

% Initialise a count
count = 1;

for tot0 = 0 : 1 : nDiagonals_vxy - 1
    for i = tot0 : -1 : 0
        j = tot0 - i;
        
        % if j1 is within the number of columns of vxy_mtrx.
        % if j2 is within the number of rows of vxy_mtrx.
        if j <= n2_k2  && i <= n1_k1
            
            % Padd the coefficients
            fxy_matrix_bi_padded_new = zero_matrix;
            fxy_matrix_bi_padded_new((i + 1):(m1)+(i + 1),(j + 1):(m2)+(j + 1)) = fxy_matrix_bi;
            
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