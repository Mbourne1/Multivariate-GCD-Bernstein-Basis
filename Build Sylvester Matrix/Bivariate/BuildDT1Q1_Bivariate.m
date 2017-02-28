function DT1Q1 = BuildDT1Q1_Bivariate(fxy, n1_k1, n2_k2)
% BuildDT1Q1_Bivariate(fxy_matrix_n,gxy_matrix_n,k1,k2,alpha,th1,th2)
%
% Build the sylvester subresultant matrix S_{k1,k2}.
%
% Inputs.
%
%
% fxy : Coefficients of the polynomial f(x,y) in the Bernstein basis
%
% n1_k1 : The degree of polynomial v(x,y) with respect to x
%
% n2_k2 : The degree of polynomial v(x,y) with respect to y
%
% % Outputs
%
% DT1Q1 : The matrix

% Get the degree of f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);

% Build two Cauchy matrices, the first for coefficients of fxy and the
% second for the coefficients of gxy
T1 = BuildT1_Bivariate(fxy, n1_k1, n2_k2);

% Build the diagonal matrix D^{-1}
D = BuildD_Bivariate_2Polys(m1, m2, n1_k1, n2_k2);

% Build the diagonal matrix Q_{1}
Q1 = BuildQ1_Bivariate(n1_k1, n2_k2);

% Build the matrix DTQ
DT1Q1 = D*T1*Q1;

end
