function [fxy_noisy, noise_matrix] = AddVariableNoiseToPoly(fxy, el, eu)
%
% Add noise to the coefficients of polynomial f(x,y)
%
% Inputs
%
% fxy : (Matrix) Coefficients of exact polynomial f(x,y)
%
% el : signal to noise low limit
%
% eu : signal to noise upper limit
%
% % Outputs
%
% fxy_noisy : (Matrix) Coefficients of the noisy form of polynomial f(x,y)
%
% noise_matrix : (Matrix) Contains noise added to the coefficients of f(x,y)


global SETTINGS
rng(SETTINGS.SEED)

% get the degree of input polynomial f
[m1, m2] = GetDegree_Bivariate(fxy);


% Get random variables

% Get set of random variables in the interval [-1, 1]
r = (2*rand(m1 + 1, m2 + 1)) - ones(m1 + 1, m2 + 1);


epslon = eu .* ones(m1 + 1, m2 + 1) -  r.*(eu - el);

% Get noise matrix
noise_matrix = r .* fxy .* epslon;

% Get f(x,y) with noise
fxy_noisy = fxy + noise_matrix;


end