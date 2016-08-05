function fy = root_y(r)
% Given a simple root r, get the vector of the polynomial f(y) whose root 
% is at r, where f(y) is a polynomial in Bernstein form.
%
%  
% Inputs.
% 
% r : root
%
% Outpus.
%
% fx :   B_{0}^{1}(y) B_{1}^{1}(y) 
%        ____________ ____________
%       | ___________|____________|
%




fy = [-r 1-r];

end