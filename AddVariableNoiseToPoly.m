function [fxy_noisy, noise_matrix] = AddVariableNoiseToPoly(fxy, el, eu)
%
% Add noise to the coefficients of polynomial f(x,y)
%
% Inputs
%
% fxy : Coefficients of exact polynomial f(x,y)
%
% el : signal to noise low limit
%
% eu : signal to noise upper limit
%
% % Outputs
%
% fxy_noisy : Coefficients of the noisy form of polynomial f(x,y)
%
% noise_matrix : Matrix of noise added to the coefficients of f(x,y)


global SETTINGS
rng(SETTINGS.SEED)

% get the degree of input polynomial f
[m1, m2] = GetDegree_Bivariate(fxy);


% Get random variables

y = (2*rand(m1+1,m2+1))-ones(m1+1,m2+1);

s = eu *ones(m1+1,m2+1) -  y.*(eu-el);

% Get noise matrix
noise_matrix = fxy.*s;

% Get f(x,y) with noise
fxy_noisy = fxy + noise_matrix;


end