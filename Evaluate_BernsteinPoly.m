
function [sum] = Evaluate_BernsteinPoly(t1,t2,fxy)

% Evaluate the bivariate bernstein basis polynomial at values x = t1, and y
% = t2

% Get the degree of the polynomial with repsect to x
m1 = size(fxy,1)-1;

% Get the degree of the polynomial with respect to y
m2 = size(fxy,2)-1;

% Perform a summing function to evaluate the curve at the given point.
sum = 0;
for i = 0:1:m1
    
    for j = 0:1:m2
        temp_val = fxy(i+1,j+1) * nchoosek(m1,i) * nchoosek(m2,j) ...
            * (1-t1)^(m1-i) * t1^i ...
            * (1-t2)^(m2-j) * t2^j;
        sum = sum + temp_val;
        
    end
    
end
end
