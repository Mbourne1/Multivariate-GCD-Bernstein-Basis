function [] = o_roots_mymethod(fxy_matrix,M)

% Get the factors of f(x,y)
[wx,vDegt_wx] = o_roots_mymethod_x(fxy_matrix,M);

LineBreakLarge()
LineBreakLarge()

[wy,vDegt_wy] = o_roots_mymethod_y(fxy_matrix,M);


[wx,wy,wxy] = o_roots_mymethod_xy(wx,wy,vDegt_wx,vDegt_wy)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

