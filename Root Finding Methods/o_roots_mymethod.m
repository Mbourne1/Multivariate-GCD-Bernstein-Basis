function [] = o_roots_mymethod(fxy, m)
% o_roots_mymethod(fxy_matrix,M)
%
% Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
% 
% m : (Int) Total degree of f(x,y) 


% Get the factors of f(x,y) with respect to x
[arr_wxy_1, vDegree_x_wxy, vDegree_y_wxy] = o_roots_mymethod_x(fxy, m);

LineBreakLarge()
LineBreakLarge()

% Get the factors of f(x,y) with respect to y
[arr_wxy_2,vDegree_y_wxy, vDegree_x_wxy] = o_roots_mymethod_y(fxy, m);

LineBreakLarge()
LineBreakLarge()

[arr_wx,arr_wy,arr_wxy] = o_roots_mymethod_xy(arr_wxy_1, arr_wxy_2)


end

