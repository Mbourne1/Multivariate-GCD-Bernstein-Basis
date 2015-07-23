function [alpha, theta1,theta2] = OptimalAlphaTheta(max_mtrx_f,min_mtrx_f, max_mtrx_g,min_mtrx_g)


% define vector f
f = [1 -1 0 0 0]

% get the degree of polynomial f and g
m1 = size(max_mtrx_f,1) -1;
m2 = size(max_mtrx_f,2) -1;
n1 = size(max_mtrx_g,1) -1;
n2 = size(max_mtrx_g,2) -1;

% Assemble the four submatrices of Matrix A
PartOne = zeros((m1+1)*(m2+1),5);
count = 1;
for i1 = 0:1:m1
    for i2 = 0:1:m2
        % Replace all of the 'count-th row'
        PartOne(count,:) = [1 0 -i1 -i2 0];
        count = count + 1 ;
    end
end

PartTwo = zeros((n1+1)*(n2+1),5);
count = 1;
for i1 = 0:1:n1
    for i2 = 0:1:n2
        % Replace all of the 'count-th row'
        PartTwo(count,:) = [1 0 -i1 -i2 -1];
        count = count + 1;
    end
end

PartThree = zeros((m1+1)*(m2+1),5);
count = 1;
for i1 = 0:1:m1
    for i2 = 0:1:m2
        % Replace all of the 'count-th row'
        PartThree(count,:) = [0 -1 i1 i2 0];
        count = count + 1;
    end
end

PartFour = zeros((n1+1)*(n2+1),5);
count = 1;
for i1 = 0:1:n1
    for i2 = 0:1:n2
        % Replace all of the 'count-th row'
        PartFour(count,:) = [0 -1 i1 i2 -1];
        count = count + 1;
    end
end

A =-[PartOne; PartTwo; PartThree; PartFour]


% Now build the vector b

lambda_vec = zeros((m1+1)*(m2+1),1);
count = 1;
for i1 = 0:1:m1
    for i2 = 0:1:m2
        lambda_vec(count) = max_mtrx_f(i1+1,i2+1);
        count = count + 1;
    end
end

mu_vec = zeros((n1+1)*(n2+1),1);
count = 1;
for i1 = 0:1:n1
    for i2 = 0:1:n2
        mu_vec(count) = max_mtrx_g(i1+1,i2+1);
        count = count + 1;
    end
end

rho_vec = zeros((m1+1)*(m2+1),1);
count = 1;
for i1 = 0:1:m1
    for i2 = 0:1:m2
        rho_vec(count) = min_mtrx_f(i1+1,i2+1);
        count = count + 1;
    end
end

tau_vec = zeros((n1+1)*(n2+1),1);
count = 1;
for i1 = 0:1:n1
    for i2 = 0:1:n2
        tau_vec(count) = min_mtrx_g(i1+1,i2+1);
        count = count + 1;
    end
end

b = -[log10(abs(lambda_vec)); log10(abs(mu_vec)); -log10(abs(rho_vec));-log10(abs(tau_vec))]


x = linprog(f,A,b)

alpha  = 10^x(5);
theta1 = 10^x(3);
theta2 = 10^x(4);



end