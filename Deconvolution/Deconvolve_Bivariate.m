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


% Get degrees of polynomials f(x,y) and g(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Build the matrix C_{m1-n1,m2-n2}(g) 
T1 = BuildT1_Bivariate(gxy, m1 - n1, m2 - n2);

% Get the polynomial f in vector form
vec_fxy = GetAsVector(fxy);

% Build the matrix D^{-1}
D = BuildD_Bivariate_2Polys(0, 0, m1, m2);

% Build the Matrix G
G = BuildG_Bivariate(m1 - n1, m2 - n2);

% Get coefficients vector h containing coefficients of polynomial h(x,y)
h = SolveAx_b(D*T1*G, vec_fxy);

% Get vector h as a matrix of coefficients of h(x,y)
hxy = GetAsMatrix(h, m1 - n1, m2 - n2);

end