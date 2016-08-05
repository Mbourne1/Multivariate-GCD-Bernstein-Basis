function fx = root_x(r)
% Given a simple root r, get the vector of the polynomial f(x) whose root 
% is at r, where f(x) is a polynomial in Bernstein form.
%
% Input.
% 
% r : root
%
% Output.
%                     ___
% fx :  B_{0}^{1}(x) |___|
%       B_{1}^{1}(x) |___|
%


%
fx = [-r; 1-r];

end