function [px] = Bernstein_Multiply(fx,gx)
% Given the Coefficients of two polynomials f(x) and g(x) in vector form, 
% output the coefficients of the product p(x).
%
% % Input.
%
% fx : Column vector of coefficients of the Bernstein Polynomial f(x)
%
% gx : Column vector of coefficients of the Bernstein Polynomial g(x)
%
% % Output.
%
% px : Column vector of coefficients of the Bernstein Polynomial p(x)

% Get the degree of polynomial f(x)
[r,~] = size(fx);
m = r - 1;

% Get the degree of polynomial g(x)
[r,~] = size(gx);
n = r - 1;

% Binomial coefficients corresponding to f(x)
Bi_m = zeros(m+1,1);
for i = 0:1:m
    Bi_m(i+1) = nchoosek(m,i);
end

Bi_n = zeros(n+1,1);
for i = 0:1:n
    Bi_n(i+1) = nchoosek(n,i);
end

fw = fx .* Bi_m;

%% Build matrix C
for i = 0:1:n
    C(i+1:(m+1)+i,i+1) = fw;
end

% Build matrix D
D_diag = zeros(1,m+n+1);
for i = 0:1:m+n
    D_diag(i+1) = 1./ nchoosek(m+n,i);
end
D = diag(D_diag);

% Build Q
Q = diag(Bi_n);

px = D*C*Q*gx;
end