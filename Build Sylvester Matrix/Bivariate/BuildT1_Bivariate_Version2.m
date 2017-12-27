function T1 = BuildT1_Bivariate_Version2(fxy, n1_k1, n2_k2)
% % Build the matrix T_{n1 - k1, n2 - k2}(f(x,y)), which is the
% (n_{1}-k_{1}, n_{2} - k_{2}) th convolution matrix of the polynomial
% f(x, y)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% n1_k1 : (Int) Degree of v(x,y) with respect to x
%
% n2_k2 : (Int) Degree of v(x,y) with respect to y
%
% % Outputs
%
% T1 : (Matrix) Convolution matrix T_{n1-k1,n2-k2}(f(x,y))


% Get degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get (m2 + 1) vectors from the matrix of coefficients f(x,y) so
% that each vector contains coefficients of a polynomail f_{i}(x)
% for i = 0,\dots,m_{2}
arrVectors = cell(m2 + 1 , 1);

for i = 0 : 1 : m2
    arrVectors{i + 1} = fxy(:, i + 1);
end

%
arrMatrix = cell(m2 + n2_k2 + 1, n2_k2 + 1);

% For each column-partition of T_{n1 - k1, n2 - k2}(f(x,y))
for j = 0 : 1 : n2_k2
    
    % For each row-partition of T_{n1 - k1, n2 - k2}(f(x,y))
    for i = 0 : 1 : m2 + n2_k2
        
        if (i - j >= 0) && (i - j <= m2)
            
            % Get f_{i}(x)
            fx = arrVectors{i - j + 1};
            
            % Build matrix T_{n1-k1}(f_{i-j}(x))
            arrMatrix{i + 1, j + 1} = BuildT1_Univariate(fx, n1_k1) .* nchoosek( m2, i - j );
            
            
        else
            
            arrMatrix{i + 1, j + 1} = zeros(m1 + n1_k1 + 1, n1_k1 + 1);
            
        end
        
        
        
        
    end
end

T1 = cell2mat(arrMatrix);

end