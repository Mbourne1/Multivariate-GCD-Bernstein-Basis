function fxy = GetWithoutThetas(fww, th1, th2)
% GetWithoutThetas(fww_matrix,th1,th2)
%
% Remove thetas from a matrix f(\omega_{1},\omega_{2}) to obtain
% coefficients of polynomial f(x,y)
%
% Inputs.
%
% fww : (Matrix) Coefficients of polynomial f(\omega_{1},\omega_{2})
%
% th1 : (Float) theta_{1}
%
% th2 : (Float) theta_{2}
%
% % Outputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)

% Get the degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fww);

% Divide the rows by theta_{1}^i and columns by theta_{2}^{j}
pre_mat  = diag(1./(th1.^(0:1:m1)));
post_mat = diag(1./(th2.^(0:1:m2)));

% 
fxy = pre_mat * fww * post_mat;

end