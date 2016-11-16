function DTQ = BuildDTQ(fxy, gxy, k1, k2)
% BuildSubresultant(fxy_matrix_n,gxy_matrix_n,k1,k2,alpha,th1,th2)
%
% Build the sylvester subresultant matrix S_{k1,k2}.
%
% % Inputs.
%
%
% fxy : Coefficients of the polynomial f(x,y)
%
% gxy : Coefficients of the polynomial g(x,y)
%
% k1 : The degree k_{1} with respect to x of the polynomial d_{k_{1},k_{2}}
%
% k2 : The degree k_{2} with respect to y of the polynomial d_{k_{1},k_{2}}
%
% % Outputs
% 
% DTQ : Matrix DTQ


% Get the degree of f(x,y) with respect to x and y
[m1,m2] = GetDegree(fxy);

% Get the degree of g(x,y) with respect to x and y
[n1,n2] = GetDegree(gxy);

% Build the matrix D_{m1+n1-k1,m2+n2-k2} * T_{n1-k1,n2-k2}(f) * Q1_{n1-k1,n2-k2}
DT1Q1 = BuildDT1Q1(fxy,n1-k1,n2-k2);

% Build the matrix D_{m1+n1-k1,m2+n2-k2} * T_{m1-k1,m2-k2}(g) * Q1_{m1-k1,m2-k2}
DT2Q2 = BuildDT1Q1(gxy,m1-k1,m2-k2);

% Build the matrix DTQ
DTQ = [DT1Q1 DT2Q2];


end
