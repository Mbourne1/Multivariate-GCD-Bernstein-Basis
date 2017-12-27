
function DT1Q1 = BuildDT1Q1_Bivariate_Version2(fxy, n1_k1, n2_k2)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% n1_k1 : (Int) Degree of v(x,y) with respect to x, and determines number of
% columns in DT(f)Q
%
% n2_k2 : (Int) Degree of v(x,y) with respect to x, and determines number of
% columns in DT(f)Q

% Get degree of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get the m2 + 1 vectors
arrVectors = cell(m2+1, 1);
for i = 0:1:m2
    arrVectors{i+1} = fxy(:, i+1);
end

%
arrMatrix = cell(m2 + n2_k2 + 1, n2_k2 + 1);

% For each column-partition of the subresultant matrix
for j = 0 : 1 : n2_k2
    
    % For each row-partition of the subresultant matrix
    for i = 0 : 1 : m2 + n2_k2

        if (i - j >= 0) && (i - j <= m2)
            
            % Get f_{i}(x)
            fx = arrVectors{i - j + 1};
            
            % Build matrix 
            arrMatrix{i + 1, j + 1} = ...
                BuildDT1Q1_Univariate(fx, n1_k1) ...
                * nchoosek(m2, i-j) ...
                * nchoosek(n2_k2, j) ...
                ./ nchoosek(m2 + n2_k2, i);
            
            
        else
            
            arrMatrix{i+1,j+1} = zeros(m1 + n1_k1 + 1, n1_k1 + 1);
            
        end
        
        
        
        
    end
end

DT1Q1 = cell2mat(arrMatrix);

end