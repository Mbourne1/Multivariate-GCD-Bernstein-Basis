function [fxy_matrix] = BuildPoly_NonSeparable(roots_f)
% Given a set of roots for f build the polynomial so that we obtain a set
% of coefficients for the polynomial in the Bernstein basis.
% 
% Inputs.
%
% roots_f : Array of factors of f(x,y)
%
% Outputs.
%
% fxy_matrix : Coefficients of the polynomial f(x,y)


% get the number of roots
[nEntries,~] = size(roots_f);
num_roots = nEntries;

% Initialise the product
prod = roots_f{1,1};

% For every other root 2,...,
for i = 2:1:num_roots
   % Multiply product by the new root
   prod = Multiply_Bivariate(prod,roots_f{i,1}) ;
end

% Output the product of all roots : f(x,y) 
fxy_matrix = prod;


end