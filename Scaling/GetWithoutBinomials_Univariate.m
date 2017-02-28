function fx = GetWithoutBinomials_Univariate(fx_bi)
% Given the polynomial coefficients f(x) in the scaled Bernstein basis, 
% divide the coefficients a_{i}\nchoosek(m,i) by nchoosek(m,i), to obtain 
% coefficients in the standard Bernstein form.

% Get the degree of polynomial f(x).
m = GetDegree_Univariate(fx_bi);

% Get the binomial coefficients corresponding to f(x)
bi_m = GetBinomials(m);

% Multiply coefficients of f(x) by the binomial coefficients
fx = fx_bi./bi_m;


end