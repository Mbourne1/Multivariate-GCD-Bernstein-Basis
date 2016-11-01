function hxy_matrix = Bernstein_Multiply_Bivariate(fxy_matrix,gxy_matrix)
% Given two Bernstein Polynomials where the coefficients are given in a 
% matrix 
% rows relate to the basis element wrt x 
% cols relate to the basis element wrt y

%               B_{0}(y)  B_{1}(y) ....
%               ________ ________
%    B_{0}(x)  |________|________| ...
%    B_{1}(x)  |________|________| ...
%    B_{2}(x)  |________|________| ...
%       ...       ...       ...

% Get the dimensions and degrees of polynomial fxy
[m1,m2] = GetDegree(fxy_matrix);

% Get the dimensions and degrees of polynomial gxy
[n1,n2] = GetDegree(gxy_matrix);

% Get the dimensions and degrees of the product hxy = fxy * gxy
o1 = m1 + n1;
o2 = m2 + n2;

% % Build the mulitplication matrix D*T*Q

DT1Q1 = BuildDT1Q1(fxy_matrix,n1,n2);

% % Build the vector of the coefficients of g(x,y)
g_vec = GetAsVector(gxy_matrix);

% % Perform multiplication to obtain hxy in vector form
h_vec = DT1Q1 * g_vec;

% % Convert hxy to matrix form
hxy_matrix = GetAsMatrix(h_vec,o1,o2);

end



