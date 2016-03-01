function [bi_m] = GetBinomials(m)
% Calculate the binomial coefficients mC0,...,mCm and return as a 
% column vector

% Initialise the column vector of binomials
bi_m = zeros(m,1);

% for each entry of the vector of binomials calculate nchoosek
for i = 0:1:m
    bi_m(i+1) = nchoosek(m,i);
end

end