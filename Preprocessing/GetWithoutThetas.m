function fxy_matrix = GetWithoutThetas(fww_matrix, th1, th2)
% GetWithoutThetas(fww_matrix,th1,th2)
%
% Remove thetas from a matrix f(\omega_{1},\omega_{2}) to obtain
% coefficients of polynomial f(x,y)
%
% Inputs.
%
% fww_matrix : Coefficients of polynomial f(\omega_{1},\omega_{2})
%
% th1 : theta_{1}
%
% th2 : theta_{2}

% Get the degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fww_matrix);

% Divide the rows by theta_{1}^i and columns by theta_{2}^{j}
pre_mat  = diag(1./(th1.^(0:1:m1)));
post_mat = diag(1./(th2.^(0:1:m2)));


fxy_matrix = pre_mat * fww_matrix * post_mat;

end