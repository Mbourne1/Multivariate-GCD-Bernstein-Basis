function DTQ = BuildDTQ_Bivariate_3Polys(fxy, gxy, hxy, k1, k2)
% BuildDTQ(fxy, gxy, hxy, k1, k2)
%
% Build the sylvester subresultant matrix S_{k1,k2}(f,g,h).
%
% % Inputs.
%
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
% 
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of the polynomial h(x,y)
%
% k1 : (Int) The degree k_{1} with respect to x of the polynomial d_{k_{1},k_{2}}
%
% k2 : (Int) The degree k_{2} with respect to y of the polynomial d_{k_{1},k_{2}}
%
% % Outputs
% 
% DTQ : (Matrix) Bivariate Sylvester subresultant matrix


% Get the degree of f(x,y), g(x,y) and h(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);

% Build the matrix D_{m1+n1-k1,m2+n2-k2} * T_{n1-k1,n2-k2}(f) * Q1_{n1-k1,n2-k2}
DT1Q1 = BuildDT1Q1_Bivariate(fxy, n1 - k1, n2 - k2);
DT2Q2 = BuildDT1Q1_Bivariate(fxy, o1 - k1, o2 - k2);

% Build the matrix D_{m1+n1-k1,m2+n2-k2} * T_{m1-k1,m2-k2}(g) * Q1_{m1-k1,m2-k2}
DT3Q3 = BuildDT1Q1_Bivariate(gxy, m1 - k1, m2 - k2);
DT4Q4 = BuildDT1Q1_Bivariate(hxy, m1 - k1, m2 - k2);

% Build the matrix DTQ
diagonal = blkdiag(DT1Q1, DT2Q2);
column = [DT3Q3; DT4Q4];

DTQ = [diagonal column];


end
