function [implic_matrix,symbolic_expression] = Implicitize_Bezier_Sylvester(CP)
% Given a set of Bezier control points for the Bezier curve. Implicitize 
% the Curve. This code corresponds to the implicitization notes in
% Sederberg CAGD notes Section 16.6.
%
% % Inputs.
%
% CP : Vector of control points stored as follows
%        ___________________
%       |x0  x1  x2 ...  xn |
%       |y0  y1  y2 ...  yn |
%       |___________________|
%
% % Outputs.
%
% implic : The coefficients of the polynomial in implicit form.

% Get the degree of the polynomial
n = size(CP,2)-1;

% get the coefficients of polynomial f
x = sym('x');

xt = CP(1,:)';
f = (xt - x);


% Multiply f by corresponding binomial coefficient
binoms_f = zeros(n+1,1);
for i = 0:1:n
    binoms_f(i+1) = nchoosek(n,i);
end

f = f .* binoms_f;

% get the coefficients of polynomial g
y = sym('y');

yt = CP(2,:)';
g = (yt - y);

g = g .* binoms_f;

% Build a zero matrix C1
C1 = zeros(2*n,n);
C1 = sym(C1);

% Build the matrix C1
for i = 1:1:n
    C1(i:i+n,i) = f;
end

% Build a zero matrix C2
C2 = zeros(2*n,n);
C2 = sym(C2);

% Build the matrix C2
for i = 1:1:n
    C2(i:i+n,i) = g;
end


% Build the matrix C2
S = [C1 C2];

%fprintf('The implicit representation is given as');
symbolic_expression = expand(det(S));





% Split the poly into a series of polynomials in terms of y
cx_vec = coeffs(symbolic_expression, x);

% a degree 3 polynomial has a bivariate implicit representation of degree
% (3,3) = (t1,t2)

max_x = n;
% get highest power of y
max_y = n;
mat = zeros(max_x+1,max_y+1);

% for each column in cx_vec
[~,c] = size(cx_vec);

for i = 1:1:c
    
    % Get the coefficients of the polynomial in y, corresponding to x^{i}
    % x^{i} * (y + y^{2} + ... )
    
    coef_vec_y = fliplr(sym2poly(cx_vec(1,i)));
    
    [~,cy] = size(coef_vec_y);
    
    mat(i,1:cy) = coef_vec_y;
    
    
end

implic_matrix = mat;

end
