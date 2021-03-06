function [fx_Bb] = PowerToBernstein_Univariate(fx)
% PowerToBernstein(fx)
%
% Given a column vector of coefficients in the power basis, convert to bernstein
% basis
% Algorithms for Polynomials in Bernstein Form
%
% Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x) in the power basis.
%
% Outputs.
%
% fx_Bb : (Vector) Coefficients of polynomail f(x) in the Bernstein basis. 


% Get the degree of f(x)
m = GetDegree_Univariate(fx);

% Initialise a vector to store coefficients in Bernstein form
fx_Bb = zeros(m+1, 1);

for j = 0:1:m
    temp_sum = 0;

    for k = 0:1:j
        temp_sum = temp_sum + ...
            (...
                nchoosek(j, k)...
                ./ nchoosek(m, k) ...
                .* fx(k+1)...
            );
    end
    fx_Bb(j+1) = temp_sum;
end


end