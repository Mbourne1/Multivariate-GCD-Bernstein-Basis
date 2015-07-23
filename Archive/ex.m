function [  ] = ex(  )
%EX Summary of this function goes here
%   Detailed explanation goes here


x = sym('x');
y = sym('y');



input = [-225 -330 -120 0 ; -330 -484 -176 0 ; -120 -176 -64 0]

[a,b] = size(input)
rows = a;
m1 = rows-1;
cols = b;
m2 = b-1;

% for each row
for i = 0:1:m1
    % for each column
    for j = 0:1:m2
        BINOMS(i+1,j+1) = nchoosek(m1,i) * nchoosek(m2,j);
    end
end

CoefficientMatrix = input./BINOMS





% Given a coefficient matrix
...CoefficientMatrix = [1 2 3; 4 5 6 ; 7 8 9; 10 11 12]
    [m2,m1] = size(CoefficientMatrix);

% let m1 be the degree in terms of x
m1 = m1-1;
m2 = m2-1;

% get max of m1 and m2
M = max(m1,m2)

% get the degree by which m1 is elevated
r = M-m1;

% get the degree by which m2 is elevated
s = M-m2;

r = 2
s = 3

% let m2 be the degree in terms of y
myPoly = 0;

for j = 0:1:m1
    for i = 0:1:m2
        myPoly = myPoly + ...
            CoefficientMatrix(i+1,j+1) ...
            *x^j * (1-x)^(m1-j) ...
            *y^i * (1-y)^(m2-i) ...
            * nchoosek(m1,j) ...
            * nchoosek(m2,i);
    end
end

expand(myPoly)

CoefficientMatrix_Full = bivariate_degree_elevation(CoefficientMatrix,r,s)'


[m2,m1] = size(CoefficientMatrix_Full);
% let m1 be the degree in terms of x
m1 = m1-1;
m2 = m2-1;
% let m2 be the degree in terms of y
myPoly = 0;

% for each column j
for j = 0:1:m1
    % for each row i
    for i = 0:1:m2
        myPoly = myPoly + ...
            CoefficientMatrix_Full(i+1,j+1) ...
            *x^j * (1-x)^(m1-j) ...
            *y^i * (1-y)^(m2-i) ...
            * nchoosek(m1,j) ...
            * nchoosek(m2,i);
    end
end
expand(myPoly)
end

