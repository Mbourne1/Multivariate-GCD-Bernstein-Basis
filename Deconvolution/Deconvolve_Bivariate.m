function [hxy_matrix] = Deconvolve_Bivariate(fxy_matrix,gxy_matrix)
% Deconvolve_Bivariate(fxy_matrix,gxy_matrix)
%
% Perform deconvolution f/g by the matrix multiplication C(g) * h = f.
%
% Inputs.
%
% fxy_matrix : Coefficients of polynomial f(x,y)
%
% gxy_matrix : Coefficients of polynomial g(x,y)
%
% Outputs.
%
% hxy_matrix : Coefficients of the polynomial h(x,y)


% Get degrees of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get the degrees of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% Build the matrix C_{m1-n1,m2-n2}(g) 
C1 = BuildT1(gxy_matrix,m1-n1,m2-n2);

% Get the polynomial f in vector form
f = GetAsVector(fxy_matrix);

% Build the matrix D^{-1}
D = BuildD(0,0,m1,m2);

% Build the Matrix G
G = BuildG(m1-n1,m2-n2);

% Get coefficients vector h containing coefficients of polynomial h(x,y)
h = SolveAx_b(D*C1*G,f);

% Get vector h as a matrix of coefficients of h(x,y)
hxy_matrix = GetAsMatrix(h,m1-n1,m2-n2);

end