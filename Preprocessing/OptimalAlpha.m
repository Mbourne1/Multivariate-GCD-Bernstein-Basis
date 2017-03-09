function [alpha] = OptimalAlpha(max_mtrx_fxy, min_mtrx_fxy, max_mtrx_gxy, min_mtrx_gxy)
% Obtain the optimal value of alpha in preprocessing the Sylvester matrix
% S_{k_{1},k_{2}}(f,g)
%
% Inputs.
%
% max_mtrx_fxy : (Matrix) Containing entries of maximum magnitude of each
% coefficient of f(x,y) in the Sylvester matrix S_{k1,k2}
%
% min_mtrx_fxy : (Matrix) Containing entries of maximum magnitude of each
% coefficient of f(x,y) in the Sylvester matrix S_{k1,k2}
%
% max_mtrx_gxy : (Matrix) Containing entries of maximum magnitude of each
% coefficient of g(x,y) in the Sylvester matrix S_{k1,k2}
%
% min_mtrx_gxy : (Matrix) Containing entries of maximum magnitude of each
% coefficient of g(x,y) in the Sylvester matrix S_{k1,k2}
%
% Outputs
%
% alpha : (Float) Optimal value \alpha such that the two partitions of 
% S_{k_{1},k_{2}} are balanced.


% Define vector f (unrelated to coefficients of polynomial f)
f = [1 -1 0];

% get the degree of polynomial f(x,y) and g(x,y)
[m1, m2] = GetDegree_Bivariate(max_mtrx_fxy);
[n1, n2] = GetDegree_Bivariate(max_mtrx_gxy);

nEntries_fxy = (m1 + 1) * (m2 + 1);
nEntries_gxy = (n1 + 1) * (n2 + 1);


% % 
% %
% %
% Assemble the four submatrices of Matrix A
PartOne = ...
    [
    ones(nEntries_fxy,1) ...
    zeros(nEntries_fxy,1) ...
    zeros(nEntries_fxy,1) ...
    ];
    

PartTwo = ...
    [
    ones(nEntries_gxy,1) ...
    zeros(nEntries_gxy,1) ...
    -1.* ones(nEntries_gxy,1)
    ];

PartThree = ...
    [
    zeros(nEntries_fxy,1) ...
    -1.* ones(nEntries_fxy,1) ...
    zeros(nEntries_fxy,1) ...
    ];

PartFour = ...
    [
    zeros(nEntries_gxy,1) ...
    -1.* ones(nEntries_gxy,1) ...
    ones(nEntries_gxy,1) ...
    ];

% Now build the vector b
lambda_vec = GetAsVector(max_mtrx_fxy);
mu_vec = GetAsVector(max_mtrx_gxy);
rho_vec = GetAsVector(min_mtrx_fxy);
tau_vec = GetAsVector(min_mtrx_gxy);


% % Find any zeros in the lambda vector
% indeces = find(~lambda_vec);
% PartOne(indeces,:) = [];
% lambda_vec(indeces,:) = [];
% 
% % Find any zeros in the mu vector
% indeces = find(~mu_vec);
% PartTwo(indeces,:) = [];
% mu_vec(indeces,:) = [];
% 
% % Find any zeros in the rho vector
% indeces = find(~rho_vec);
% PartThree(indeces,:) = [];
% rho_vec(indeces,:) = [];
% 
% % Find any zeros in the tau vector
% indeces = find(~tau_vec);
% PartFour(indeces,:) = [];
% tau_vec(indeces,:) = [];


b = [log10(lambda_vec); log10(mu_vec); -log10(rho_vec);-log10(tau_vec)]';

A = [PartOne; PartTwo; PartThree; PartFour];



x = linprog(f,-A,-b);

try
    alpha  = 10^x(3);
catch
    alpha = 1;
    return
end


end