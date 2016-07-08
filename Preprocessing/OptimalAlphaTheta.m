function [alpha, theta1,theta2] = OptimalAlphaTheta(max_mtrx_f,min_mtrx_f, max_mtrx_g,min_mtrx_g)
% Obtain the optimal value of alpha and theta

% define vector f
f = [1 -1 0 0 0];

% Get the degree of polynomial f(x,y)
[m1,m2] = GetDegree(max_mtrx_f);

% Get the degree of polynomial g(x,y)
[n1,n2] = GetDegree(max_mtrx_g);

% Assemble the four submatrices of Matrix A

nEntries_f = (m1 + 1) * (m2 + 1);
nEntries_g = (n1 + 1) * (n2 + 1);


v_i1_f = GetAsVector(diag(0:1:m1) * ones(m1+1,m2+1));
v_i2_f = GetAsVector(ones(m1+1,m2+1) * diag(0:1:m2));

v_i1_g = GetAsVector(diag(0:1:n1) * ones(n1+1,n2+1));
v_i2_g = GetAsVector(ones(n1+1,n2+1) * diag(0:1:n2));

PartOne = ...
    [
    ones(nEntries_f,1) ...
    zeros(nEntries_f,1) ...
    -1 .* v_i1_f ...
    -1 .* v_i2_f ...
    zeros(nEntries_f,1)
    ];

PartTwo = ...
    [
    ones(nEntries_g,1) ...
    zeros(nEntries_g,1) ...
    -1.* v_i1_g ...
    -1.* v_i2_g ...
    -1.* ones(nEntries_g,1)
    ];

PartThree = ...
    [
    zeros(nEntries_f,1) ...
    -1.* ones(nEntries_f,1) ...
    v_i1_f ...
    v_i2_f ...
    zeros(nEntries_f,1)
    ];

PartFour = ...
    [
    zeros(nEntries_g,1) ...
    -1 .* ones(nEntries_g,1) ...
    v_i1_g ...
    v_i2_g ...
    ones(nEntries_g,1)
    ];

% Now build the vector b
lambda_vec = GetAsVector(max_mtrx_f);
mu_vec = GetAsVector(max_mtrx_g);
rho_vec = GetAsVector(min_mtrx_f);
tau_vec = GetAsVector(min_mtrx_g);

% % Find any zeros in the lambda vector
indeces = find(~lambda_vec);
PartOne(indeces,:) = [];
lambda_vec(indeces,:) = [];

% Find any zeros in the mu vector
indeces = find(~mu_vec);
PartTwo(indeces,:) = [];
mu_vec(indeces,:) = [];

% Find any zeros in the rho vector
indeces = find(~rho_vec);
PartThree(indeces,:) = [];
rho_vec(indeces,:) = [];

% Find any zeros in the tau vector
indeces = find(~tau_vec);
PartFour(indeces,:) = [];
tau_vec(indeces,:) = [];


b = [log10(lambda_vec); log10(mu_vec); -log10(rho_vec);-log10(tau_vec)]';

A = [PartOne; PartTwo; PartThree; PartFour];

x = linprog(f,-A,-b);

try
    alpha  = 10^x(5);
    theta1 = 10^x(3);
    theta2 = 10^x(4);
catch
    alpha = 1;
    theta1 = 1;
    theta2 = 1;


end

end