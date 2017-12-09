function Q = BuildQ_Bivariate_2Polys(n1_k1, n2_k2, m1_k1, m2_k2)
% BuildQ_Bivariate(n1_k1, n2_k2, m1_k1, m2_k2)
%
% Build the diagonal matrix Q whose entries contain binomial coefficients
% corresponding to the entries of u(x,y) and v(x,y) in the vector [v;-u]
% such that Q*[v;-u] gives the coefficients of u(x,y) and v(x,y) in scaled
% Bernstein form. 
% The matrix Q forms a part of the modified kth Sylvester subresultant
% matrix D^{-1}[T(f) T(g)]Q.
%
% Inputs.
%
% n1_k1 : (Int) Degree of polynomial v(x,y) with respect to x
%
% n2_k2 : (Int) Degree of polynomial v(x,y) with respect to y
%
% m1_k1 : (Int) Degree of polynomial u(x,y) with respect to x
%
% m2_k2 : (Int) Degree of polynomial u(x,y) with respect to y

% Build the matrix Q_{1}
Q1 = BuildQ1_Bivariate(n1_k1, n2_k2);

% Build the matrix Q_{2}
Q2 = BuildQ1_Bivariate(m1_k1, m2_k2);

% Build the matrix Q
Q = blkdiag(Q1, Q2);

end