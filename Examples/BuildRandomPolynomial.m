function [f_roots] = BuildRandomPolynomial(m, interval_low, interval_high)
%
% Build random univariate polynomial by producing a matrix of roots and
% corresponding multiplicities.
%
% % Inputs
%
% m : (Int) Degree of polynomial f(x)
%
% interval_low : (Float) : Upper bound for all roots
%
% interval_up : (Float) : Lower bound for all roots
%
% % Outputs
%
% f_roots : (Matrix) Two column matrix of the form [ root | multiplicity]


global SETTINGS

a = interval_low;
b = interval_high;

% Get a multiplicity structure for polynomial f
prob_arr = zeros(1,m);
for i = 1:1:m
    prob_arr(i) = i./ nchoosek(m+1,2);
end

prob_arr = fliplr(prob_arr);

rng(SETTINGS.SEED);

% Get the multiplicity structure of d.
total = 0;
i = 1;
while total < m
    r = rand;
    prob = prob_arr;
    x = sum(r >= cumsum([0, prob]));
    if (total + x) <= m
        mult_arr_d(i) = x;
        total = total + x;
        i = i+1;
    end
end

% get the number of roots of d
nRoots_t = length(mult_arr_d);


% Get the roots
% Get a set of unique roots
% the 1000 and 1000 contain the roots to the unit interval
detail = 100;
format 'long';



roots = a + randperm(detail,nRoots_t)./(detail./(b-a));

f_roots = [roots' mult_arr_d'];
end