function [m1, m2] = GetDegree_Bivariate(fxy)
% GetDegree(fxy)
%
% Get the degree strucuture m_{1} and m_{2} of the bivariate polynomial
% f(x,y)
%
% Inputs.
%
% fxy : Matrix of coefficients of polynomial f(x,y)
%

% Get the dimensions of the matrix of coefficients of f(x,y)
[r,c] = size(fxy);

% Get the degree with respect to x.
m1 = r - 1;

% Get the degree with respect to y.
m2 = c - 1;

end