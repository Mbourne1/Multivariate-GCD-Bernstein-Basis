function G = BuildG_Bivariate(t1, t2)
% BuildG_Bivariate(t1, t2)
%
% Build the matrix G used in the APF 
%
% % Inputs
%
% t1 : Degree of d(x,y) with respect to x
%
% t2 : Degree of d(x,y) with respect to y 

% H [C(u); C(v)] G d = [f;g]



G = BuildQ1_Bivariate(t1,t2);

end