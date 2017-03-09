function [DTQ] = BuildDTQ_Bivariate_2Polys_Version2(fxy, gxy, k1, k2)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% k1 : (Int)
%
% k2 : (Int)


% Get degree of f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of g(x,y) with respect to x and y
[n1, n2] = GetDegree_Bivariate(gxy);

% Build matrices D^{-1}T(f)Q and D^{-1}T(g)Q
DT1Q1 = BuildDT1Q1_Bivariate_2Polys_Version2(fxy, n1-k1, n2-k2);
DT2Q2 = BuildDT1Q1_Bivariate_2Polys_Version2(gxy, m1-k1, m2-k2);

% Build DTQ
DTQ = [DT1Q1 DT2Q2];



end


function DT1Q1 = BuildDT1Q1_Bivariate_2Polys_Version2(fxy, n1_k1, n2_k2)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% n1_k1 : Degree of v(x,y) with respect to x, and determines number of
% columns in DT(f)Q
%
% n2_k2 : Degree of v(x,y) with respect to x, and determines number of
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

for j = 0:1:n2_k2
    
    for i = 0:1:m2+n2_k2

        if (i-j>=0) && (i-j <= m2)
            
            % Get f_{i}(x)
            fx = arrVectors{i-j+1};
            
            % Build matrix 
            arrMatrix{i+1,j+1} = BuildDT1Q1_Univariate(fx, n1_k1) * nchoosek(m2, i-j) * nchoosek(n2_k2, j) ./ nchoosek(m2+n2_k2, i);
            
            
        else
            
            arrMatrix{i+1,j+1} = zeros(m1+n1_k1+1, n1_k1+1);
            
        end
        
        
        
        
    end
end

DT1Q1 = cell2mat(arrMatrix);

end