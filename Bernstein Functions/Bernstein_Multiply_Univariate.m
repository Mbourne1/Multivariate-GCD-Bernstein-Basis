function [px] = Bernstein_Multiply_Univariate(fx, gx)
% Given the Coefficients of two polynomials f(x) and g(x) in vector form, 
% output the coefficients of the product p(x).
%
% % Input.
%
% fx : (Vector) Vector of coefficients of the Bernstein Polynomial f(x)
%
% gx : (Vector) Vector of coefficients of the Bernstein Polynomial g(x)
%
% % Output.
%
% px : (Vector) Column vector of coefficients of the Bernstein Polynomial p(x)

if( size(fx, 2) > 1 || size(gx, 2) > 1)
   error('fx and gx must be column vectors'); 
end


% Get the degree of polynomial f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Binomial coefficients corresponding to f(x)
Bi_m = GetBinomials(m);

% Get scaled 
fx_bi = fx .* Bi_m;

% % Build matrix C
T_fx = zeros(m + n + 1, n + 1);

for i = 0:1:n
    T_fx(i + 1 : (m + 1) + i, i + 1) = fx_bi;
end

D = BuildD_Univariate_2Polys(m, n);
Q = BuildQ1_Univariate(n);

% Get the product of f(x) and g(x)
px = D * T_fx * Q * gx;


end