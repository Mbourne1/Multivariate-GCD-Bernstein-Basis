function [DTQ] = BuildDTQ_Univariate(fx, gx, k)
% BuildDTQ_Univariate(fx, gx, t)
%
% Build the matrix DTQ = D^{-1}T(f,g)*Q.
%
% Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x) 
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% k : (Int) Index of subresultant
%
% Note: If you wish to build the Sylvester matrix for preprocessed and
% scaled polynomials f(\omega)g(\omega), the preprocessed forms must be the
% inputs to this function.


% Get degree of polynomials f(x) and g(x)
m = GetDegree_Univariate(fx);
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