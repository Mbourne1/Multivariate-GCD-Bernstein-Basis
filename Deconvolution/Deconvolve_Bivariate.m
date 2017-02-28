function [hxy] = Deconvolve_Bivariate(fxy, gxy)
% Deconvolve_Bivariate(fxy_matrix,gxy_matrix)
%
% Perform deconvolution f/g by the matrix multiplication C(g) * h = f.
%
% Inputs.
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% Outputs.
%
% hxy : (Matrix) Coefficients of the polynomial h(x,y)


% Get degrees of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get the degrees of polynomial g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Build the matrix C_{m1-n1,m2-n2}(g) 
C1 = BuildT1_bivar(gxy, m1-n1, m2-n2);

% Get the polynomial f in vector form
vec_f = GetAsVector(fxy);

% Build the matrix D^{-1}
D = BuildD_bivar(0,0,m1,m2);

% Build the Matrix G
G = BuildG_bivar(m1-n1,m2-n2);

% Get coefficients vector h containing coefficients of polynomial h(x,y)
h = SolveAx_b(D*C1*G,vec_f);

% Get vector h as a matrix of coefficients of h(x,y)
hxy = GetAsMatrix(h,m1-n1,m2-n2);

end