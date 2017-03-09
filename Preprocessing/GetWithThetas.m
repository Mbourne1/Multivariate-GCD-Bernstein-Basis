function fww = GetWithThetas(fxy, th1, th2)
% Given the matrix in f(x,y) obtain the matrix of f(\theta_{1},\theta_{2})
%
% Inputs.
%
% fxy : (Matrix) Coefficients of f(x,y)
% 
% th1 : Optimal theta_{1}
% 
% th2 : Optimal theta_{2}
%
% % Outputs
%
% fww : (Matrix) Coefficients of f(\omega_{1}, \omega_{2}0

% Get the degree structure of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get f(w,w) by multiplying rows by \theta_{1} and cols by \theta_{2}.
fww = diag(th1.^(0:1:m1)) * fxy * diag(th2.^(0:1:m2));



end