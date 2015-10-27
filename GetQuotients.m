function [uxy_matrix_calc, vxy_matrix_calc,...
    opt_theta_1, opt_theta_2] = GetQuotients(fxy_matrix_working, gxy_matrix_working,...
    t1,t2)
% Initialise Global Variables
global bool_Q
global bool_preproc

%%
% Get the degrees of polynomial f(x,y)
[r,c] = size(fxy_matrix_working);
m1 = r - 1;
m2 = c - 1;

% Get the degrees of polynomial g(x,y)
[r,c] = size(gxy_matrix_working);
n1 = r - 1;
n2 = c - 1;


%%
fprintf('Now Build the corresponding subresultant\n')

switch bool_preproc
    case 'y'
        
        [fxy_max_mtrx,fxy_min_mtrx] = GetMaxMin(fxy_matrix_working,n1,n2,t1,t2);
        [gxy_max_mtrx,gxy_min_mtrx] = GetMaxMin(gxy_matrix_working,m1,m2,t1,t2);
        
        
        [opt_theta_1, opt_theta_2] = ...
            OptimalTheta(fxy_max_mtrx,fxy_min_mtrx,gxy_max_mtrx,gxy_min_mtrx);
    case 'n'
        opt_theta_1 = 1;
        opt_theta_2 = 1;
end
% Build the t-th subresultant
St = BuildSubresultant(fxy_matrix_working,gxy_matrix_working,t1,t2,opt_theta_1,opt_theta_2);


%% Find Optimal column for removal from St
% given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized

[~,cols] = size(St);

% QR Decomposition of the Sylvester Matrix S_{k}
[Qk,Rk] = qr(St);
n = n1+n2;
for j=1:1:cols
    
    ck = St(:,j);
    [Q,~] = qrdelete(Qk,Rk,j);
    cd = Q'*ck;
    d = cd(n+1:end,:);
    residuals_QR(j) = norm(d);
    
end

%Obtain the column for which the residual is minimal.
[~,opt_col] = min(log10(residuals_QR));
fprintf('Optimal column for removal is given by %i \n',opt_col)


%% Get the coefficients for u(x,y) and v(x,y)

Atj = St;
cki = St(:,opt_col);
Atj(:,opt_col) = [];

[~,n_col] = size(Atj);
[Q,R] = qr(Atj);
R1 = R(1:n_col,:);
cd = Q'*cki;
c = cd(1:n_col,:);
x_ls = R1\c;

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ]

num_coeff_v = (n1-t1+1) * (n2-t2+1);
num_coeff_u = (m1-t1+1) * (m2-t2+1);


switch bool_Q
    case 'y'
        % Q is included in the Sylvester matrix
        % uw and vw exclude binomial coefficients
        
        % get coefficients of u and v
        vw_calc = vecx(1:num_coeff_v);
        uw_calc = -vecx(num_coeff_v+1:end);
        
    case 'n'
        % Q is excluded from the Sylvester matrix, and included in the
        % vector x
        % uw and vw include binomial coefficients.
        
        % get coefficients of u and v
        vw_calc = vecx(1:num_coeff_v);
        uw_calc = -vecx(num_coeff_v+1:end);
        
        %%
        % Remove binomial coefficients from vw and uw matrices
        [r,c] = size(vw_calc)
        n1_t1 = r - 1;
        n2_t2 = c - 1;
        
        bi_n1_t1 = zeros(n1_t1,1);
        for i = 0:1:n1_t1
            bi_n1_t1(i+1) = nchoosek(n1_t1,i);
        end
        bi_n2_t2 = zeros(n2_t2,1);
        for i = 0:1:n2_t2
            bi_n2_t2(i+1) = nchoosek(n2_t2,i);
        end
        
        % divide the rows of vw by nchoosek(n1_t1,i)
        vw_calc = diag(1./bi_n1_t1) * vw_calc
        % divide the cols of vw by nchoosek(n2_t2,j)
        vw_calc = vw_calc * diag(1./bi_n2_t2)
        
        %%
        % Remove the binomial coefficients from uw
        [r,c] = size(uw_calc)
        m1_t1 = r - 1;
        m2_t2 = c - 1;
        
        bi_m1_t1 = zeros(m1_t1,1);
        for i = 0:1:m1_t1
            bi_m1_t1(i+1) = nchoosek(m1_t1,i);
        end
        bi_m2_t2 = zeros(m2_t2,1);
        for i = 0:1:m2_t2
            bi_m2_t2(i+1) = nchoosek(m2_t2,i);
        end
        
        % divide the rows of vw by nchoosek(n1_t1,i)
        uw_calc = diag(1./bi_m1_t1) * uw_calc;
        % divide the cols of uw by nchoosek(n2_t2,j)
        uw_calc = uw_calc * diag(1./bi_m2_t2);
        
        
        
end


%% Obtain u(x,y) in its matrix form
% Arrange uw into a matrix form based on its dimensions.

% Initialise a counter
count = 0;

% Get number of diagaonals in u(x,y) matrix
num_diags = max(m1-t1+m2-t2);

% for each anti-diagonal
for tot = 0:1:num_diags
    for i = tot:-1:0
        j = tot - i;
        if i <= m1-t1 && j <= m2-t2
            uw_calc_mtrx(i+1,j+1) = uw_calc(count+1);
            count = count + 1;
        end
    end
end

%% Obtain v(x,y) in its matrix form
% Arrange vw into a matrix form based on their dimensions.
count = 0;

% Get number of diagonals in v(x,y) matrix
num_diags = max(n1-t1+n2-t2);

% For each anti-diagonal
for tot = 0:1:num_diags
    for i = tot:-1:0
        j = tot - i;
        if i <= n1-t1 && j <= n2-t2
            vw_calc_mtrx(i+1,j+1) = vw_calc(count+1);
            count = count + 1;
        end
    end
end

%%
% Remove the thetas from the matrix of coefficients of v(x,y)

% Initialise an empty matrix to store vx - the coefficients of polynomial v
% with thetas removed.
vxy_matrix_calc = zeros(n1-t1+1,n2-t2+1);

% Remove the theta1 from the rows for each row
for i1 = 0:1:(n1-t1)
    vxy_matrix_calc(i1+1,:) = vw_calc_mtrx(i1+1,:) ./ (opt_theta_1^i1);
end
% Remove the theta2 from the columns
for i2 = 0:1:(n2-t2)
    vxy_matrix_calc(:,i2+1) = vxy_matrix_calc(:,i2+1) ./ (opt_theta_2^i2);
end

%%
% Remove the thetas from the matrix of coefficients of u(x,y)

% for each row, divide by theta2^i1
uxy_matrix_calc = zeros(m1-t1+1,m2-t2+1);

% Remove the theta1 from the rows for each row
for i1 = 0:1:m1-t1
    uxy_matrix_calc(i1+1,:) = uw_calc_mtrx(i1+1,:) ./ (opt_theta_1^i1);
end
% Remove the theta2 from the columns
for i2 = 0:1:m2-t2
    uxy_matrix_calc(:,i2+1) = uxy_matrix_calc(:,i2+1) ./ (opt_theta_2^i2);
end


end