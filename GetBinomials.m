function [bi_m] = GetBinomials(m)
% Calculate the binomial coefficients mC0,...,mCm and return as a 
% column vector
%
% % Inputs
%
% m : (Int) Degree of polynomial 
%
% % Outputs
%
% bi_m : Binomials (mC0, mC1, ..., mCm)

% Initialise the column vector of binomials
bi_m = zeros(m,1);

% For each entry of the vector of binomials calculate nchoosek
for i = 0:1:m
    
    bi_m(i+1,1) = nchoosek(m,i);
    
end

end