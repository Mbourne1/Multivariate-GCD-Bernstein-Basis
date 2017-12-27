function hxy = Bernstein_Multiply_Bivariate(fxy, gxy)
% Bernstein_Multiply_Bivariate(fxy_matrix, gxy_matrix)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% % Outputs
%
% hxy : (Matrix) Coefficients of polynomial h(x,y), the product of f(x,y)
% and g(x,y)
%

% Get the dimensions and degrees of polynomial f(x, y) and g(x, y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Get the dimensions and degrees of the product hxy = fxy * gxy
o1 = m1 + n1;
o2 = m2 + n2;

% Build the mulitplication matrix D*T*Q
DT1Q1 = BuildDT1Q1_Bivariate(fxy, n1, n2);

% Build the vector of the coefficients of g(x,y)
v_gxy = GetAsVector(gxy);

% Perform multiplication to obtain hxy in vector form
h_vec = DT1Q1 * v_gxy;

% Convert hxy to matrix form
hxy = GetAsMatrix(h_vec, o1, o2);

end



