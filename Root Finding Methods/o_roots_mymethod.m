function [] = o_roots_mymethod(fxy_matrix,m)
% o_roots_mymethod(fxy_matrix,M)
%
% Inputs
%
% fxy_matrix : Coefficients of polynomial f(x,y)
% 
% m : degree of f(x,y)


% Get the factors of f(x,y) with respect to x
[wx,vDegt_wx] = o_roots_mymethod_x(fxy_matrix,m);

LineBreakLarge()
LineBreakLarge()

% Get the factors of f(x,y) with respect to y
[wy,vDegt_wy] = o_roots_mymethod_y(fxy_matrix,m);

LineBreakLarge()
LineBreakLarge()

[wx,wy,wxy] = o_roots_mymethod_xy(wx,wy,vDegt_wx,vDegt_wy)


end

