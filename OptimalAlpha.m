function [alpha] = OptimalAlpha(max_mtrx_f,min_mtrx_f, max_mtrx_g,min_mtrx_g)
% Obtain the optimal value of alpha

% define vector f
f = [1 -1 0];

% get the degree of polynomial f and g

[m1,m2] = GetDegree(max_mtrx_f);
[n1,n2] = GetDegree(max_mtrx_g);

nEntries_f = (m1 + 1) * (m2 + 1);
nEntries_g = (n1 + 1) * (n2 + 1);


% % 
% %
% %
% Assemble the four submatrices of Matrix A
PartOne = ...
    [
    ones(nEntries_f,1) ...
    zeros(nEntries_f,1) ...
    zeros(nEntries_f,1) ...
    ];
    

PartTwo = ...
    [
    ones(nEntries_g,1) ...
    zeros(nEntries_g,1) ...
    -1.* ones(nEntries_g,1)
    ];

PartThree = ...
    [
    zeros(nEntries_f,1) ...
    -1.* ones(nEntries_f,1) ...
    zeros(nEntries_f,1) ...
    ];

PartFour = ...
    [
    zeros(nEntries_g,1) ...
    -1.* ones(nEntries_g,1) ...
    ones(nEntries_g,1) ...
    ];

% Now build the vector b
lambda_vec = GetAsVector(max_mtrx_f);
mu_vec = GetAsVector(max_mtrx_g);
rho_vec = GetAsVector(min_mtrx_f);
tau_vec = GetAsVector(min_mtrx_g);


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