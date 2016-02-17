function mat = root_x(r)
% Given a simple root r, get the matrix of its coefficients 
%
%                 ___
%   B_{0}^{1}(x) |___|
%   B_{1}^{1}(x) |___|
%
%
% Input.
% 
%   r : root
%

%%
mat = [-r; 1-r];

end