function [hxy] = Bernstein_Multiply_Bivariate(fxy,gxy)
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
[m1,m2] = GetDegree(fxy);

% Get the dimensions and degrees of polynomial gxy
[n1,n2] = GetDegree(gxy);

DTQ = BuildDT1Q1(fxy,n1,n2);
g = GetAsVector(gxy);

h = DTQ * g;

hxy = GetAsMatrix(h,m1+n1,m2+n2);

end