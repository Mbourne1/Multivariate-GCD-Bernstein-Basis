function [alpha, th1, th2] = OptimalAlphaTheta(max_mtrx_f,min_mtrx_f, max_mtrx_g,min_mtrx_g)
% OptimalAlphaTheta(max_mtrx_f,min_mtrx_f, max_mtrx_g,min_mtrx_g)
%
% Obtain the optimal value of alpha, \theta_{1} and \theta_{2} in
% preprocessing the Sylvester matrix
% S_{k_{1},k_{2}}(f,g)
%
% Inputs.
%
% max_mtrx_f : Matrix containing entries of maximum magnitude of each
% coefficient of f(x,y) in the Sylvester matrix S_{k1,k2}
%
% min_mtrx_f : Matrix containing entries of maximum magnitude of each
% coefficient of f(x,y) in the Sylvester matrix S_{k1,k2}
%
% max_mtrx_g : Matrix containing entries of maximum magnitude of each
% coefficient of g(x,y) in the Sylvester matrix S_{k1,k2}
%
% min_mtrx_g : Matrix containing entries of maximum magnitude of each
% coefficient of g(x,y) in the Sylvester matrix S_{k1,k2}
%
% Outputs
%
% alpha : Optimal value \alpha such that the two partitions of
% S_{k_{1},k_{2}} are balanced.
%
% th1 : Optiml value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}

% Define vector f
f = [1 -1 0 0 0];

% Get the degree of polynomial f(x,y)
[m1,m2] = GetDegree_Bivariate(max_mtrx_f);

% Get the degree of polynomial g(x,y)
[n1,n2] = GetDegree_Bivariate(max_mtrx_g);

% Assemble the four submatrices of Matrix A

nCoeffs_f = (m1 + 1) * (m2 + 1);
nCoeffs_g = (n1 + 1) * (n2 + 1);


v_i1_f = GetAsVector(diag(0:1:m1) * ones(m1+1,m2+1));
v_i2_f = GetAsVector(ones(m1+1,m2+1) * diag(0:1:m2));

v_i1_g = GetAsVector(diag(0:1:n1) * ones(n1+1,n2+1));
v_i2_g = GetAsVector(ones(n1+1,n2+1) * diag(0:1:n2));

PartOne = ...
    [
    ones(nCoeffs_f,1) ...
    zeros(nCoeffs_f,1) ...
    -1 .* v_i1_f ...
    -1 .* v_i2_f ...
    zeros(nCoeffs_f,1)
    ];

PartTwo = ...
    [
    ones(nCoeffs_g,1) ...
    zeros(nCoeffs_g,1) ...
    -1.* v_i1_g ...
    -1.* v_i2_g ...
    -1.* ones(nCoeffs_g,1)
    ];

PartThree = ...
    [
    zeros(nCoeffs_f,1) ...
    -1.* ones(nCoeffs_f,1) ...
    v_i1_f ...
    v_i2_f ...
    zeros(nCoeffs_f,1)
    ];

PartFour = ...
    [
    zeros(nCoeffs_g,1) ...
    -1 .* ones(nCoeffs_g,1) ...
    v_i1_g ...
    v_i2_g ...
    ones(nCoeffs_g,1)
    ];

% Now build the vector b
lambda_vec = GetAsVector(abs(max_mtrx_f));
mu_vec = GetAsVector(abs(max_mtrx_g));
rho_vec = GetAsVector(abs(min_mtrx_f));
tau_vec = GetAsVector(abs(min_mtrx_g));

% % Find any zeros in the lambda vector
indeces = find(lambda_vec==0);
PartOne(indeces,:) = [];
lambda_vec(indeces,:) = [];

% Find any zeros in the mu vector
indeces = find(mu_vec==0);
PartTwo(indeces,:) = [];
mu_vec(indeces,:) = [];

% Find any zeros in the rho vector
indeces = find(rho_vec==0);
PartThree(indeces,:) = [];
rho_vec(indeces,:) = [];

% Find any zeros in the tau vector
indeces = find(tau_vec==0);
PartFour(indeces,:) = [];
tau_vec(indeces,:) = [];


b = [log10(lambda_vec); log10(mu_vec); -log10(rho_vec);-log10(tau_vec)];

A = [PartOne; PartTwo; PartThree; PartFour];

warning('off')
x = linprog(f,-A,-b);
warning('on')
try
    
    th1 = 10^x(3);
    th2 = 10^x(4);
    alpha  = 10^x(5);
catch
    alpha = 1;
    th1 = 1;
    th2 = 1;
end

end