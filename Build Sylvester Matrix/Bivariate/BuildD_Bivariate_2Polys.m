function D = BuildD_Bivariate_2Polys(m1, m2, n1_k1, n2_k2)
% BuildD_Bivariate_2Polys(m1, m2, n1_k1, n2_k2)
%
% Build the matrix D
%
% % Inputs
%
% m1 : (Int) Degree of polynomial f(x,y) with respect to x
%
% m2 : (Int) Degree of polynomial f(x,y) with respect to y
%
% n1_k1 : (Int) Degree of polynomial v(x,y) with respect to x
%
% n2_k2 : (Int) Degree of polynomial v(x,y) with respect to y
%
% % Outputs
%
% D : (Matrix) Matrix D^{-1}_{m1+n1-k1,m2+n2-k2}

%
D_mat = GetWithBinomials_Bivariate(ones(m1 + n1_k1 + 1, m2 + n2_k2 + 1));


% Get the D_matrix as a vector
D_vec = GetAsVector_Version1(1./D_mat);


% Form a diagonal matrix from the vector D_vec
D = diag(D_vec);

end

