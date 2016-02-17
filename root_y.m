function mat = root_y(r)
% Given a simple root r, get the matrix of its coefficients of the factor 
% (x-r) in Bernstein form : -rB_{0} + (1-r)B_{1} 
%
%                
%   B_{0}^{1}(y) B_{1}^{1}(y) 
%   _________________________
%  | ___________|____________|
%
%
% Inputs.
% 
% r : root
%
%%

mat = [-r 1-r];

end