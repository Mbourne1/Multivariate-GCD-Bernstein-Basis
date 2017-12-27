function [] = DifferentiatePowerBasis(f)
% Given a set of coefficients in the Power basis (Where the leading 
% coefficient is of highest degree. Differentiate to obtain a set of
% coefficients of the derivative.


% Get the degree of the polynomial.
m = size(f,1) -1;

% multiply each coefficient by the old power
temp_vec = m:-1:0;
temp_vec = temp_vec';

f_der = temp_vec .* f
f_der(m+1) = []
% Remove the last coefficient



end