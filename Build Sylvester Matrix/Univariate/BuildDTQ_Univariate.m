function [DTQ] = BuildDTQ_Univariate(fx, gx, k)
% BuildDTQ_Univariate(fx, gx, t)
%
% Build the matrix DTQ = D^{-1}T(f,g)*Q.
%
% Inputs.
%
% [fx, gx] : (Vector) Coefficients of polynomial f(x) and g(x)
%
% k : index of subresultant
%
% Note: If you wish to build the Sylvester matrix for preprocessed and
% scaled polynomials f(\omega)g(\omega), the preprocessed forms must be the
% inputs to this function.


% Get degree of polynomial f(w)
m = GetDegree_Univariate(fx);

% Get degree of polynomial g(w)
n = GetDegree_Univariate(gx);

% Build matrix D^{-1}
D = BuildD_Univariate_2Polys(m, n-k);

% Build matrix T(f,g) = T_{n-k}(f) T_{m-k}*(g)
T = BuildT_Univariate_2Polys(fx, gx, k);

% Build matrix Q = [Q_{n-k} Q_{m-k}]
Q = BuildQ_Univariate_2Polys(m, n, k);

% Get D^{-1} * T(f,g) * Q
DTQ = D*T*Q;



end