function [fxy_matrix] = BuildPoly_NonSeparable(roots_f)
% given a set of roots for f build the polynomial so that we obtain a set
% of coefficients for the polynomial in the Bernstein basis.

% get the number of roots
[r,c] = size(roots_f);
num_roots = r;

prod = roots_f{1,1};

for i = 2:1:num_roots
   prod = BernMultiply_Bivariate(prod,roots_f{i,1}) ;
end
fxy_matrix = prod;
end