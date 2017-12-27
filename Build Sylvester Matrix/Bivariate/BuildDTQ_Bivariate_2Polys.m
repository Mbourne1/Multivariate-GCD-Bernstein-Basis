function DTQ = BuildDTQ_Bivariate_2Polys(fxy, gxy, k1, k2)
% BuildDTQ_Bivariate_2Polys(fxy, gxy, k1, k2)
%
% Build the sylvester subresultant matrix S_{k1,k2}.
%
% % Inputs.
%
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% k1 : (Int) The degree k_{1} with respect to x of the polynomial d_{k_{1},k_{2}}
%
% k2 : (Int) The degree k_{2} with respect to y of the polynomial d_{k_{1},k_{2}}
%
% % Outputs
% 
% DTQ : (Matrix) D^{-1}T(f(x,y),g(x,y))Q


% Get the degree of f(x,y) and g(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Build the matrix D_{m1+n1-k1,m2+n2-k2} * T_{n1-k1,n2-k2}(f) * Q1_{n1-k1,n2-k2}
DT1Q1 = BuildDT1Q1_Bivariate(fxy, n1 - k1, n2 - k2);

% Build the matrix D_{m1+n1-k1,m2+n2-k2} * T_{m1-k1,m2-k2}(g) * Q1_{m1-k1,m2-k2}
DT2Q2 = BuildDT1Q1_Bivariate(gxy, m1 - k1, m2 - k2);

% Build the matrix DTQ
DTQ = [DT1Q1 DT2Q2];


end
