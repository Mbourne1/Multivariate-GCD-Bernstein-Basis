function [theta1,theta2] = OptimalTheta(max_matrix_fxy, min_matrix_fxy, max_matrix_gxy, min_matrix_gxy)
% OptimalTheta(max_mtrx_f,min_mtrx_f, max_mtrx_g,min_mtrx_g)
%
% Obtain the optimal value of \theta_{1} and \theta_{2} in
% preprocessing the Sylvester matrix
% S_{k_{1},k_{2}}(f,g)
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
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}

% Define vector f (Note : Not associated with polynomial f)
f = [1 -1 0 0];

% Get the degree of polynomial f(x,y) and g(x,y)
[m1, m2] = GetDegree_Bivariate(max_matrix_fxy);
[n1, n2] = GetDegree_Bivariate(max_matrix_gxy);

% Assemble the four submatrices of Matrix A
PartOne = zeros((m1+1)*(m2+1),4);
count = 1;
for i1 = 0:1:m1
    for i2 = 0:1:m2
        % Replace all of the 'count-th row'
        PartOne(count,:) = [1 0 -i1 -i2 ];
        count = count + 1 ;
    end
end

PartTwo = zeros((n1+1)*(n2+1),4);
count = 1;
for i1 = 0:1:n1
    for i2 = 0:1:n2
        % Replace all of the 'count-th row'
        PartTwo(count,:) = [1 0 -i1 -i2 ];
        count = count + 1;
    end
end

PartThree = zeros((m1+1)*(m2+1),4);
count = 1;
for i1 = 0:1:m1
    for i2 = 0:1:m2
        % Replace all of the 'count-th row'
        PartThree(count,:) = [0 -1 i1 i2 ];
        count = count + 1;
    end
end

PartFour = zeros((n1+1)*(n2+1),4);
count = 1;
for i1 = 0:1:n1
    for i2 = 0:1:n2
        % Replace all of the 'count-th row'
        PartFour(count,:) = [0 -1 i1 i2 ];
        count = count + 1;
    end
end

A =[PartOne; PartTwo; PartThree; PartFour];


% Now build the vector b

lambda_vec = zeros((m1+1)*(m2+1),1);
count = 1;
for i1 = 0:1:m1
    for i2 = 0:1:m2
        lambda_vec(count) = max_matrix_fxy(i1+1,i2+1);
        count = count + 1;
    end
end

mu_vec = zeros((n1+1)*(n2+1),1);
count = 1;
for i1 = 0:1:n1
    for i2 = 0:1:n2
        mu_vec(count) = max_matrix_gxy(i1+1,i2+1);
        count = count + 1;
    end
end

rho_vec = zeros((m1+1)*(m2+1),1);
count = 1;
for i1 = 0:1:m1
    for i2 = 0:1:m2
        rho_vec(count) = min_matrix_fxy(i1+1,i2+1);
        count = count + 1;
    end
end

tau_vec = zeros((n1+1)*(n2+1),1);
count = 1;
for i1 = 0:1:n1
    for i2 = 0:1:n2
        tau_vec(count) = min_matrix_gxy(i1+1,i2+1);
        count = count + 1;
    end
end

b = [log10(abs(lambda_vec)); log10(abs(mu_vec)); -log10(abs(rho_vec));-log10(abs(tau_vec))];


x = linprog(f,-A,-b);

try
    theta1 = 10^x(3);
    theta2 = 10^x(4);
catch
    theta1 = 1;
    theta2 = 1;
end


end