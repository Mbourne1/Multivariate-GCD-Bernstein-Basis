function [AM_fx] = GetArithmeticMean(fxy, n1_k1, n2_k2)
%
% Compute the arithmetic mean of the entries of f(x,y) in the convolution
% matrix C_{n_{1}-k_{1},n_{2}-k_{2}}(f(x,y))
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% n1_k1 : (Int) 
%
% n2_k2 : (Int) 
%
% % Outputs
%
% AM_fx : (Float)


[m1, m2] = GetDegree_Bivariate(fxy);

AM_fx = ...
    (m1 + n1_k1 + 1)  ...
    * (m2 + n2_k2 + 1) ...
    / (m1+1)^2 ...
    / (m2+1)^2 ...
    / (n1_k1 + 1) ...
    / (n2_k2 + 1) ...
    * sum(fxy(:));

    
end