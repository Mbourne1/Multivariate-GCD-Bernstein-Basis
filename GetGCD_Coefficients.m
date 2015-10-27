function [dxy_calc_mtrx] = GetGCD_Coefficients(uxy_calc_matrix,vxy_calc_matrix,...
    fxy_matrix_working,gxy_matrix_working,...
    t1,t2,...
    opt_theta1,opt_theta2)

% Get the degrees of polynomial f
[r,c] = size(fxy_matrix_working);
m1 = r - 1;
m2 = c - 1;

% Get the degrees of polynomial g
[r,c] = size(gxy_matrix_working);
n1 = r - 1;
n2 = c - 1;

%% Assemble the Subresultant S_{t1,t2}

% Build matrix H
H = BuildH(m1,m2,n1,n2);

% % Build Matrix C
% Build Matrix C1
C1_u = BuildC1(uxy_calc_matrix,t1,t2,m1,m2,opt_theta1,opt_theta2);
% Build Matrix C2
C2_v = BuildC1(vxy_calc_matrix,t1,t2,n1,n2,opt_theta1,opt_theta2);

C = [ C1_u;
    C2_v ];

% Buid matrix G
G = BuildG(t1,t2);

% Build the Coefficient Matrix HCG 
HCG = H*C*G;

%% Include thetas in f(x,y) to obtain f(\theta_1, \theta_2)
% Build Vector f
% Get f with thetas included
fxy_matrix_th = zeros(m1+1,m2+1);
for i1 = 0:1:m1
    fxy_matrix_th(i1+1,:) = fxy_matrix_working(i1+1,:) .* (opt_theta1^i1);
end
for i2 = 0:1:m2
    fxy_matrix_th(:,i2+1) = fxy_matrix_th(:,i2+1) .* (opt_theta2^i2);
end

%% Get f(x,y) as a vector

fxy_vec_th = getAsVector(fxy_matrix_th)

%% Include thetas in g(x,y) to obtain g(\theta_1, theta_2)

% Build Vector g of coefficients of polynomial g
gxy_matrix_th = zeros(n1+1,n2+1);
for i1 = 0:1:n1
    gxy_matrix_th(i1+1,:) = gxy_matrix_working(i1+1,:) .* (opt_theta1^i1);
end
for i2 = 0:1:n2
    gxy_matrix_th(:,i2+1) = gxy_matrix_th(:,i2+1) .* (opt_theta2^i2);
end

%% Get g(x,y) as a vector

gxy_vec_th = getAsVector(gxy_matrix_th)


%% Create the right hand side vector
rhs_vec = [...
    fxy_vec_th;
    gxy_vec_th];


%% Obtain vector x
x = pinv(HCG) * rhs_vec;
dw_calc = x;

residual = pinv(HCG)*rhs_vec - x

%%
% Arrange dw into a matrix form based on its dimensions.
dw_calc_mtrx = getAsMatrix(dw_calc,t1,t2)

%%
% remove thetas from dw
% Divide the row i by theta1^i
for i1 = 0:1:t1
    dxy_calc_mtrx(i1+1,:) = dw_calc_mtrx(i1+1,:) ./ (opt_theta1^i1);
end
for i2 = 0:1:t2
    dxy_calc_mtrx(:,i2+1) = dxy_calc_mtrx(:,i2+1) ./ (opt_theta2^i2);
end

