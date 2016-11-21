function brn_sym = Brn_Sym(CP)
% Given a set of control points of f(x) generate a symbolic expression for the
% Bernstein polynomial f(x). 

% Get the degree of the polynomial
[~,c] = size(CP);
m = c - 1;

% Initialise the Symbols B0, ..., Bm to represent the basis elements.
B = sym('B%d', [m+1 1]);

sum = 0;
for i = 0:1:m
   sum = sum + B(i+1) * CP(i+1) ;
end

brn_sym = sum;

end