function DT1Q1 = BuildDT1Q1(fxy_matrix,n1_k1,n2_k2)
% BuildSubresultant(fxy_matrix_n,gxy_matrix_n,k1,k2,alpha,th1,th2)
%
% Build the sylvester subresultant matrix S_{k1,k2}.
%
% Inputs.
%
%
% fxy_matrix_n : Coefficients of the polynomial f(x,y)
%
% k1 : The degree k_{1} with respect to x of the polynomial d_{k_{1},k_{2}}
%
% k2 : The degree k_{2} with respect to y of the polynomial d_{k_{1},k_{2}}
%
%

% Get the degree of f(x,y) with respect to x and y
[m1,m2] = GetDegree(fxy_matrix);

% Build two Cauchy matrices, the first for coefficients of fxy and the
% second for the coefficients of gxy
T1 = BuildT1(fxy_matrix,n1_k1,n2_k2);

% Build the diagonal matrix D^{-1}
D = BuildD(m1,m2,n1_k1,n2_k2);

% Build the diagonal matrix Q_{1}
Q1 = BuildQ1(n1_k1,n2_k2);


DT1Q1 = D*T1*Q1;

end
