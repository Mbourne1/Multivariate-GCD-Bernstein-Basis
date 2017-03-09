function fx = root_x(r)
% Given a simple root r, get the vector of the polynomial f(x) whose root 
% is at r, where f(x) is a polynomial in Bernstein form.
%
% Input.
% 
% r : root
%
% Output.
%               
% fx : (Vector) Coefficients of f(x) 


% Get vector of 
fx = [-r; 1-r];

end