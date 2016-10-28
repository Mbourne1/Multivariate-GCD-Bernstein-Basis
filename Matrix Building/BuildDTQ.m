function DTQ = BuildDTQ(fxy_matrix,gxy_matrix,k1,k2)
% BuildSubresultant(fxy_matrix_n,gxy_matrix_n,k1,k2,alpha,th1,th2)
%
% Build the sylvester subresultant matrix S_{k1,k2}.
%
% Inputs.
%
%
% fxy_matrix_n : Coefficients of the polynomial f(x,y)
%
% gxy_matrix_n : Matrix of coefficients of polynomial gxy in the scaled 
%                bernstein basis (including the binomial coefficients) b_{i,j}
%                \binom{n_{1}}{i}\binom{n_{2}}{j}
%
% k1 : The degree k_{1} with respect to x of the polynomial d_{k_{1},k_{2}}
%
% k2 : The degree k_{2} with respect to y of the polynomial d_{k_{1},k_{2}}



% Get the degree of f(x,y) with respect to x and y
[m1,m2] = GetDegree(fxy_matrix);

% Get the degree of g(x,y) with respect to x and y
[n1,n2] = GetDegree(gxy_matrix);


DT1Q1 = BuildDT1Q1(fxy_matrix,n1-k1,n2-k2);
DT2Q2 = BuildDT1Q1(gxy_matrix,m1-k1,m2-k2);

DTQ = [DT1Q1 DT2Q2];


end
