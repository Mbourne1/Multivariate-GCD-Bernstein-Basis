function [sum] = Evaluate_BernsteinPoly_Bivariate(t1, t2, fxy)
% Evaluate_BernsteinPoly_Bivariate(t1,t2,fxy)
%
% Evaluate the bivariate bernstein basis polynomial at values x = t1, and y
% = t2
%
% % Inputs
%
% t1 : (Int) Evaluation point x_{i}
%
% t2 : (Int) Evaluation point y_{i}
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% % Outputs
%
% sum : Value of function f(x,y) at (t_{1},t_{2})

% Get degree of f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);


% Perform a summing function to evaluate the curve at the given point.
sum = 0;

% For each row of coefficients in f(x,y)
for i = 0:1:m1
    % For each column of coefficients in f(x,y)
    for j = 0:1:m2
        temp_val = fxy(i+1,j+1) * nchoosek(m1,i) * nchoosek(m2,j) ...
            * (1-t1)^(m1-i) * t1^i ...
            * (1-t2)^(m2-j) * t2^j;
        sum = sum + temp_val;
        
    end
    
end
end
