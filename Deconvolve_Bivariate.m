function [hxy_matrix] = Deconvolve_Bivariate(fxy_matrix,gxy_matrix)
% Perform deconvolution f/g by the matrix multiplication C(g) * h = f.
%


% Get degrees of polynomial f
[m1,m2] = GetDegree(fxy_matrix);

% Get the degrees of polynomial g
[n1,n2] = GetDegree(gxy_matrix);

% Build the matrix C(g) 
C1 = BuildT1(gxy_matrix,m1-n1,m2-n2);

% Get the polynomial f in vector form
f = GetAsVector(fxy_matrix);

D = BuildD(0,0,m1,m2,0,0);
G = BuildG(m1-n1,m2-n2);

h = SolveAx_b(D*C1*G,f);

hxy_matrix = GetAsMatrix(h,m1-n1,m2-n2);

end