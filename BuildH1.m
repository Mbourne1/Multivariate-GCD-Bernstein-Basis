function H1 = BuildH1(m1,m2)
% Build H1, the diagonal matrix of binomials coefficients.
%
%                               Inputs.
%
% m1 :  Degree of polynomial f with respect to x
%
% m2 :  Degree of polynomial f with respect to y
%

Pre_matrix = diag(GetBinomials(m1))
Post_matrix = diag(GetBinomials(m2))

Binom_matrix =  Pre_matrix * ones(m1+1,m2+1) * Post_matrix ;
Binom_vector = GetAsVector(Binom_matrix);

H1 = diag(1./Binom_vector);

end