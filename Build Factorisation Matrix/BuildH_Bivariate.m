function H = BuildH_Bivariate(m1, m2, n1, n2)
% BuildH_Bivariate(m1, m2, n1, n2)
%
% Build the matrix H^{-1}, which is used in the APF function where
% H^{-1}[T_{u} ; T_{v}] [d] = [f;g] 
%
% % Inputs:
%
% m1 : (Int) Degree of polynomial f(x,y) with respect to x
%
% m2 : (Int) Degree of polynomial f(x,y) with respect to y
%
% n1 : (Int) Degree of polynomial g(x,y) with respect to x
%
% n2 : (Int) Degree of polynomial g(x,y) with respect to y
%
% % Outputs.
%
% H : (Matrix) The matrix H^{-1}
%

% Build the first partition H1
H1 = BuildH1_Bivariate(m1, m2);

% Build the second partition H2
H2 = BuildH1_Bivariate(n1, n2);

% Build the block diagonal matrix [H1 0 ; 0 H2] = diag([H1 H2])
H = blkdiag(H1,H2);
end