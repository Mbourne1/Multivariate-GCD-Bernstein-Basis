function T1 = BuildT1_Symbolic_Bivariate(m1, m2, n1_k1, n2_k2)
% BuildT1(fxy_matrix,n1_k1,n2_k2)
%
% Build the matrix T1, the first partition of the Sylvester matrix
% S_{k_{1},k_{2}}. T_{n_{1}-k_{1},n_{2}-k_{2}}(f) * [v] gives a vector of coefficients of
% the polynomial multiplication f(x,y) * v(x,y).
%
% % Inputs.
%
% m1 : (Int) Degree of polynomial f(x,y) with respect to x
%
% m2 : (Int) Degree of polynomial f(x,y) with respect to y
%
% n1_k1 : (Int) Degree of polynomial v(x,y) with respect to x
%
% n2_k2 : (Int) Degree of polynomial v(x,y) with repsect to y
%
% % Outputs
%
% T1 : (Matrix) T_{n1-k1,n2-k2}(f)

% Get symbolic polynomial
fxy = sym('A', [m1 + 1 m2 + 1]);

% Get the degree of f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);

% Get coefficients of f(x,y) in scaled Bernstein form
fxy_bi = GetWithBinomials_Bivariate(fxy);

% Initalise a zero matrix
zero_matrix = sym(zeros(m1 + n1_k1 + 1, m2 + n2_k2 + 1));

% Get the number of rows in T_{k_{1},k_{2}}(f)
nRowsT1 = (m1 + n1_k1 +1 ) * (m2 + n2_k2 +1);
nColsT1 = (n1_k1 + 1) * (n2_k2 + 1);

% Initialise a zero matrix T_{n_{1}-k_{1},n_{2}-k_{2}}
T1 = sym(zeros(nRowsT1, nColsT1));

% For every diagonal of the matrix vxy.
nDiagonals_vxy = (n2_k2 + 1) + (n1_k1 + 1) -1;

% Initialise a count
count = 1;

column_order = 'Rearranged';

switch column_order
    case 'Standard'
        for tot0 = 0:1:nDiagonals_vxy - 1
            for i = tot0 : -1:0
                j = tot0-i;
                % if j1 is within the number of columns of vxy_mtrx.
                % if j2 is within the number of rows of vxy_mtrx.
                if j <= n2_k2  && i <= n1_k1
                    
                    % Padd the coefficients
                    fxy_bi_padded_new = zero_matrix;
                    fxy_bi_padded_new(...
                        (i+1):(m1)+(i+1), ...
                        (j+1):(m2)+(j+1)) = fxy_bi;
                    
                    temp_vec = GetAsVector_Symbolic(fxy_bi_padded_new);
                    
                    % Assign column
                    T1(:,count) = temp_vec;
                    
                    % Increment Counter
                    count = count + 1;
                end
                %
            end
        end
        
    case 'Rearranged'
        
        
        for i = 0 : 1 : n1_k1
            for j = 0 : 1 : n2_k2
                
                % Padd the coefficients
                fxy_bi_padded_new = zero_matrix;
                
                fxy_bi_padded_new(...
                    (i + 1) : m1 + (i + 1),...
                    (j + 1) : m2 + (j + 1)) = fxy_bi;
                
                temp_vec = GetAsVector_Symbolic(fxy_bi_padded_new);
                
                % Assign column
                T1(:,count) = temp_vec;
                
                % Increment Counter
                count = count + 1;
            end
            %
        end
        
        
end


end