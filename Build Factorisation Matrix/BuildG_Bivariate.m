function G = BuildG_Bivariate(t1, t2)
% BuildG_Bivariate(t1, t2)
%
% Build the matrix G used in the APF 
%
% % Inputs
%
% t1 : (Int) Degree of d(x,y) with respect to x
%
% t2 : (Int) Degree of d(x,y) with respect to y 
%
% Outputs
%
% G : (Matrix) 

% Used in constructing the factorisation matrix H [C(u); C(v)] G 

% Build G
G = BuildQ1_Bivariate(t1, t2);

end