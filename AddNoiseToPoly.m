function [fxy_noisy, noise_matrix] = AddNoiseToPoly(fxy, el)
% Add noise to the coefficients of the polynomial f(x,y)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% el : (Float) Noise level
%
% % Outputs
%
% fxy_noisy : (Matrix) Coefficients of f(x,y) with noise added
%
% noise_matrix : (Matrix) Matrix of noise added to coefficients of f(x,y)


% Initialise global variables
global SETTINGS

% Set seed
rng(SETTINGS.SEED)


% Get the degrees of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Obtain a matrix of randomly distributed numbers between [-1 and 1]
rp = ( 2*rand(m1+1, m2+1)) - ones(m1+1, m2+1);

% multiply by the noise:signal ratio so so that we have a range of
% values between [-noise : + noise]
s = rp *el;

% Get the noise matrix
noise_matrix = fxy.*s;

% add the noise matrix to the exact matrix
fxy_noisy = fxy + noise_matrix;


end