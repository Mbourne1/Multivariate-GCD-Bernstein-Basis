function [m] = GetDegree_Univariate(fx)
% GetDegree(fx)
%
% Get the degree m  of the univariate polynomial
% f(x,y)
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% % Outputs
%
% m : Degree of polynomail f(x)

% Get the dimensions of the matrix of coefficients of f(x,y)
[r,c] = size(fx);

if c > 1 
   error([mfilename ' : Input must be a column vector. \n']); 
end

% Get the degree with respect to x.
m = r - 1;

end