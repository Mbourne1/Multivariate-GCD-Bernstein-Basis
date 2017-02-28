function H = BuildH_Bivariate_3Polys(m1, m2, n1, n2, o1, o2)
% BuildH_Bivariate(m1, m2, n1, n2)
%
% Build the matrix H^{-1}, which is used in the APF function where
% H^{-1}[T_{u} ; T_{v}] [d] = [f;g] 
%
% % Inputs:
%
% [m1, m2] - Degree of polynomial f with respect to x and y
%
% [n1, n2] - Degree of polynomial g with respect to x and y
%
% [o1, o2] - Degree of polynomial h with respect to x and y
%
% % Outputs.
%
% H : The matrix H^{-1}
%

% Build the first partition H1
H1 = BuildH1_Bivariate(m1, m2);

% Build the second partition H2
H2 = BuildH1_Bivariate(n1, n2);

% Build the third parititon H3
H3 = BuildH1_Bivariate(o1, o2);

% Build the block diagonal matrix [H1 0 ; 0 H2] = diag([H1 H2])
H = blkdiag(H1,H2,H3);
end