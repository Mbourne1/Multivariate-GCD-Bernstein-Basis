function [t, opt_theta_1, opt_theta_2] = GetDegree_Total(fxy_matrix,gxy_matrix,m,n)

global bool_preproc

% Get degrees of polynomial f(x,y)
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get degrees of polynomial g(x,y)
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;


% Degree elevate fxy
fxy_matrix = DegreeElevateToTotalDegree(fxy_matrix,m);
gxy_matrix = DegreeElevateToTotalDegree(gxy_matrix,n);

% Get degrees of polynomial f(x,y)
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get degrees of polynomial g(x,y)
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

min_sin_val_vec = [];
cond_vec = [];
opt_theta_1_vec = [];
opt_theta_2_vec = [];

% for each possible total degree
for k=1:1:min(m,n)
    
    %% Apply preprocessing
    switch bool_preproc
        case 'y'
            [max_mtrx_f, min_mtrx_f] = GetMaxMin(fxy_matrix,n1,n2,k,k);
            [max_mtrx_g, min_mtrx_g] = GetMaxMin(gxy_matrix,m1,m2,k,k);
            
            [opt_theta_1, opt_theta_2] = OptimalTheta(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g);
        case 'n'
            opt_theta_1 = 1;
            opt_theta_2 = 1;
            
    end
    
    opt_theta_1_vec(k) = opt_theta_1;
    opt_theta_2_vec(k) = opt_theta_2;
    
    
    %% Build the Sylvester matrix S_{k,k}
    C_f = BuildT1(fxy_matrix,n1,n2,k,k,opt_theta_1,opt_theta_2);
    C_g = BuildT1(gxy_matrix,m1,m2,k,k,opt_theta_1,opt_theta_2);
    
    Sk = [C_f C_g];
    
    % Get SVD of SK
    min_sin_val_vec(k) = min(svd(Sk));
    
    % Get the condition of Sk
    cond_vec(k) = cond(Sk);
    
end


% plot the minimum singular values
figure(1)
title('minimum Singular Value for each subresultant matrix S_{k,k}')
hold on
plot(log10(min_sin_val_vec),'-s');
xlabel('k : index of subresultant')
ylabel('log_{10} Minimum Singular Value')
hold off

figure(2)
title('Condition of each subresultant S_{k,k}')
hold on
plot(log10(cond_vec),'-s');
xlabel('k : index of subresultant S_{k}')
ylabel('log_{10} Condition Number')
hold off


% Get the max change in min singular values
[val, index] = max(diff(min_sin_val_vec));

% Set the degree of the GCD to be equal to the index of the largest change
t = index-1;

% Set the optimal theta 1 and theta 2
opt_theta_1 = opt_theta_1_vec(index-1);
opt_theta_2 = opt_theta_2_vec(index-1);

end