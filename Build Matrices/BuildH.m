function H = BuildH(m1,m2,n1,n2)
% BuildH(m1,m2,n1,n2)
%
% Build the matrix H^{-1}, which is used in the APF function where
% H^{-1}[T_{u} ; T_{v}] [d] = [f;g] 
%
% Inputs:
%
% m1 - Degree of polynomial f with respect to x
%
% m2 - Degree of polynomial f with respect to y
%
% n1 - Degree of polynomial g with respect to x
%
% n2 - Degree of polynomial g with respect to y
%
% Outputs.
%
% H : The matrix H^{-1}
%

% Build the first partition H1
H1 = BuildH1(m1,m2);

% Build the second partition H2
H2 = BuildH1(n1,n2);

% Build the block diagonal matrix [H1 0 ; 0 H2] = diag([H1 H2])
H = blkdiag(H1,H2);
end