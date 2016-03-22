function H1 = BuildH1(m1,m2)
% BuildH1(m1,m2)
%
% Build the diagonal matrix H1 consisting of binomials coefficients
%
%                               Inputs.
%
% m1 :  Degree of polynomial f with respect to x
%
% m2 :  Degree of polynomial f with respect to y
%

% Get Matrix of binomial coefficients
Binom_matrix = GetWithBinomials(ones(m1+1,m2+1));

% Convert to vector
Binom_vector = GetAsVector(Binom_matrix);

% Diagonalize
H1 = diag(1./Binom_vector);

end