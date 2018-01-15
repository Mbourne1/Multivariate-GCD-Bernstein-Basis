function [alpha, th1, th2] = OptimalAlphaTheta(max_matrix_fxy, ...
    min_matrix_fxy, max_matrix_gxy, min_matrix_gxy)
% OptimalAlphaTheta(max_mtrx_f,min_mtrx_f, max_mtrx_g,min_mtrx_g)
%
% Obtain the optimal value of alpha, \theta_{1} and \theta_{2} in
% preprocessing the Sylvester matrix
% S_{k_{1},k_{2}}(f(x,y),g(x,y)) 
%
% Inputs.
%
% max_mtrx_f : (Matrix) containing entries of maximum magnitude of each
% coefficient of f(x,y) in the Sylvester matrix S_{k1,k2}
%
% min_mtrx_f : (Matrix) containing entries of maximum magnitude of each
% coefficient of f(x,y) in the Sylvester matrix S_{k1,k2}
%
% max_mtrx_g : (Matrix) containing entries of maximum magnitude of each
% coefficient of g(x,y) in the Sylvester matrix S_{k1,k2}
%
% min_mtrx_g : (Matrix) containing entries of maximum magnitude of each
% coefficient of g(x,y) in the Sylvester matrix S_{k1,k2}
%
% Outputs
%
% alpha : (Float) Optimal value \alpha such that the two partitions of
% S_{k_{1},k_{2}}(f,g) are balanced.
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}

% Define vector f
f = [1 -1 0 0 0];

% Get the degree of polynomials f(x,y) and g(x,y)
[m1, m2] = GetDegree_Bivariate(max_matrix_fxy);
[n1, n2] = GetDegree_Bivariate(max_matrix_gxy);

% Assemble the four submatrices of Matrix A

nCoefficients_fxy = (m1 + 1) * (m2 + 1);
nCoefficients_gxy = (n1 + 1) * (n2 + 1);


v_i1_f = GetAsVector_Version1(diag(0:1:m1) * ones(m1 + 1, m2 + 1));
v_i2_f = GetAsVector_Version1(ones(m1 + 1,m2 + 1) * diag(0 : 1 : m2));

v_i1_g = GetAsVector_Version1(diag(0:1:n1) * ones(n1 + 1,n2 + 1));
v_i2_g = GetAsVector_Version1(ones(n1 + 1, n2 + 1) * diag(0 : 1 : n2));

matPartOne = ...
    [
    ones(nCoefficients_fxy, 1) ...
    zeros(nCoefficients_fxy, 1) ...
    -1 .* v_i1_f ...
    -1 .* v_i2_f ...
    zeros(nCoefficients_fxy, 1)
    ];

matPartTwo = ...
    [
    ones(nCoefficients_gxy, 1) ...
    zeros(nCoefficients_gxy, 1) ...
    -1.* v_i1_g ...
    -1.* v_i2_g ...
    -1.* ones(nCoefficients_gxy, 1)
    ];

matPartThree = ...
    [
    zeros(nCoefficients_fxy, 1) ...
    -1.* ones(nCoefficients_fxy, 1) ...
    v_i1_f ...
    v_i2_f ...
    zeros(nCoefficients_fxy, 1)
    ];

matPartFour = ...
    [
    zeros(nCoefficients_gxy,1) ...
    -1 .* ones(nCoefficients_gxy,1) ...
    v_i1_g ...
    v_i2_g ...
    ones(nCoefficients_gxy,1)
    ];

% Now build the vector b
vLambda = GetAsVector_Version1(abs(max_matrix_fxy));
vMu = GetAsVector_Version1(abs(max_matrix_gxy));
vRho = GetAsVector_Version1(abs(min_matrix_fxy));
vTau = GetAsVector_Version1(abs(min_matrix_gxy));

% Remove zeros from vectors and remove corresponding rows from matrices
[vLambda, matPartOne] = RemoveZeros(vLambda, matPartOne);
[vMu, matPartTwo] = RemoveZeros(vMu, matPartTwo);
[vRho, matPartThree] = RemoveZeros(vRho, matPartThree);
[vTau, matPartFour] = RemoveZeros(vTau, matPartFour);

% % Set up linprog problem

b = [log10(vLambda); log10(vMu); -log10(vRho);-log10(vTau)];

A = [matPartOne; matPartTwo; matPartThree; matPartFour];

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


function [vLambda, matPartOne] = RemoveZeros(vLambda, matPartOne)
% Remove any zeros in the vector and remove corresponding rows in the
% matrix
%
% vLambda : (Vector)
%
% matPartOne : (Matrix)


% % Find any zeros in the lambda vector
indeces = find(vLambda == 0);
matPartOne(indeces,:) = [];
vLambda(indeces,:) = [];

end