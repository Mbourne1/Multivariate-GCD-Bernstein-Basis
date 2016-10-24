function [C] = PowerToBernstein(fx)
% PowerToBernstein(fx)
%
% Given a column vector of coefficients in the power basis, convert to bernstein
% basis
% Algorithms for Polynomials in Bernstein Form
%
% Inputs.
%
% fx : Vector of coefficients of polynomial f(x) in the power basis.
%
% Outputs.
%
% c : Vector of coefficients of polynomail f(x) in the Bernstein basis. 


% Get the degree of f(x)
m = GetDegree(fx);

C = zeros(m+1,1);

for j = 0:1:m
    temp_sum = 0;
    for k = 0:1:j
        temp_sum = temp_sum + ...
            (...
                nchoosek(j,k)...
                ./ nchoosek(m,k) ...
                .* fx(k+1)...
            );
    end
    C(j+1) = temp_sum;
end


end