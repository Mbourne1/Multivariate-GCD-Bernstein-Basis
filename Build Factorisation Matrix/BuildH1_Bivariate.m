function H1 = BuildH1_Bivariate(m1, m2)
% BuildH1_Bivariate(m1, m2)
%
% Build the diagonal matrix H1 consisting of binomials coefficients,
% Similar to function BuildD(). Used in computing d in the matrix vector 
% problem H*[C(u);C(v)]*G *d = [f;g] 
%
% Inputs.
%
% m1 : (Int) Degree of polynomial f with respect to x
%
% m2 : (Int) Degree of polynomial f with respect to y
%
% Outputs.
%
% H1 : (Matrix) H_{1} is a matrix of binomial coefficients corresponding to the
% polynomial f(x,y) such that H_{1}*[f] gives the coefficients in scaled
% Bernstein form. 


% Get Matrix of binomial coefficients
Binom_matrix = GetWithBinomials_Bivariate(ones(m1 + 1, m2 + 1));

% Convert to vector
Binom_vector = GetAsVector(Binom_matrix);

% Diagonalize
H1 = diag(1 ./ Binom_vector);

end