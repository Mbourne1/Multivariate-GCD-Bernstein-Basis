function [] = o1(ex_num,BOOL_PREPROC)
% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                 Global Variables

% bool_q :  Whether the diagonal matrix Q is included in the Sylvester
%           Matrix or the solution vector x.
%           1 : Include in The Sylvester matrix
%           0 : Include in the solution vector
global bool_Q
bool_Q = 1;


global bool_GM
global bool_thetas

if (BOOL_PREPROC == 1)
    bool_GM = 1;
    bool_thetas = 1;
elseif (BOOL_PREPROC == 0)
    bool_GM = 0;
    bool_thetas = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Inputs

%   ex_num


% take a set of polynomial roots
[f_x_roots, f_y_roots, g_x_roots, g_y_roots, d_x_roots, d_y_roots,...
    u_x_roots, u_y_roots,v_x_roots,v_y_roots] = Examples(ex_num);

% Build Polynomial of the GCD in x terms and y terms
d_x_poly = BuildPoly(d_x_roots);
d_y_poly = BuildPoly(d_y_roots);

% Get degree of d in terms of x
d1 = length(d_x_poly)-1;

% Get degree of d in terms of y
d2 = length(d_y_poly)-1;

% Get Tdeg = Total Degree of d
d = d1 + d2;

% Build Polynomial f in x terms and y terms
f_x_poly = BuildPoly(f_x_roots);
f_y_poly = BuildPoly(f_y_roots);

% Get degree of f in terms of x
m1 = length(f_x_poly)-1;

% Get degree of f in terms of y
m2 = length(f_y_poly)-1;

% Get Tdeg of polynomial f
m = m1 + m2;

% Build Polynomial g in x terms and y terms
g_x_poly = BuildPoly(g_x_roots);
g_y_poly = BuildPoly(g_y_roots);

% Get degree of f in terms of x
n1 = length(g_x_poly)-1;

% Get degree of f in terms of y
n2 = length(g_y_poly)-1;

% Get Tdeg of polynomial g
n = n1+n2;

% Build Polynomial u in x terms and y terms
u_x_poly = BuildPoly(u_x_roots);
u_y_poly = BuildPoly(u_y_roots);

% Build Polynomial v in x terms and y terms
v_x_poly = BuildPoly(v_x_roots);
v_y_poly = BuildPoly(v_y_roots);

% get matrix of polynomial f coefficients including the binomial
% coefficients. The entries of fxy_matrix_bi are of the modified bernstein
% basis
fxy_matrix_bi =  f_x_poly' * f_y_poly

% for each row of fxy_matrix_bi. degree elevate

% get matrix of polynomial g coefficients including the binomial
% coefficients.
gxy_matrix_bi =  g_x_poly' * g_y_poly;

% get matrix of polynomial coefficients including the binomial coefficients
dxy_matrix_bi = d_x_poly' * d_y_poly;
uxy_matrix_bi = u_x_poly' * u_y_poly;
vxy_matrix_bi = v_x_poly' * v_y_poly;

% Strip the binomial coefficients from f
for i = 0:1:m1
    fxy_matrix(i+1,:) = fxy_matrix_bi(i+1,:) ./ nchoosek(m1,i);
end
for j = 0:1:m2
    fxy_matrix(:,j+1) = fxy_matrix(:,j+1) ./ nchoosek(m2,j);
end

fxy_matrix./fxy_matrix_bi

% Strip the binomial coefficients from g
for i = 0:1:n1
    gxy_matrix(i+1,:) = gxy_matrix_bi(i+1,:) ./ nchoosek(n1,i);
end
for j = 0:1:n2
    gxy_matrix(:,j+1) = gxy_matrix(:,j+1) ./ nchoosek(n2,j);
end

% Strip the binomial coefficients from u
for i = 0:1:m1-d1
    uxy_matrix(i+1,:) = uxy_matrix_bi(i+1,:) ./ nchoosek(m1-d1,i);
end
for j = 0:1:m2-d2
    uxy_matrix(:,j+1) = uxy_matrix(:,j+1) ./ nchoosek(m2-d2,j);
end

% Strip the binomial coefficients from v
for i = 0:1:n1-d1
    vxy_matrix(i+1,:) = vxy_matrix_bi(i+1,:) ./ nchoosek(n1-d1,i);
end
for j = 0:1:n2-d2
    vxy_matrix(:,j+1) = vxy_matrix(:,j+1) ./ nchoosek(n2-d2,j);
end

% TODO - add noise to the coefficients of fxy stored in matrix fxy_matrix

% TODO - add noise to the coefficients of gxy stored in matrix gxy_matrix



% given the entries fxy_matrix_bi, get the polynomial in the power basis
% f_exp = getInPowerBasis(fxy_matrix_bi);
% g_exp = getInPowerBasis(gxy_matrix_bi);
% d_exp = getInPowerBasis(dxy_matrix_bi);
% u_exp = getInPowerBasis(uxy_matrix_bi);
% v_exp = getInPowerBasis(vxy_matrix_bi);
%
% fxy_matrix_bi;
% gxy_matrix_bi;
% dxy_matrix_bi;
% uxy_matrix_bi;
% vxy_matrix_bi;

fprintf('Exact Degree of GCD : %i \n', d)
fprintf('Exact t1 : %i \n',d1)
fprintf('Exact t2 : %i \n',d2)





sing_val_vec = zeros(min(m,n),1);
t_val_vec = zeros(min(m,n),2);





for k = 1:1:min(m,n)
    min_sing_val_vec = [];
    temp_vec_2 = [];
    count = 1;
    
    for k1 = 0:1:k
        k2 = k-k1;
        if (k1 > m1 || k1 > n1) || (k2 > m2 || k2 > n2)
        else
            
            % STAGE 1 - Preprocess by geometric mean
      
            switch bool_GM
                case 1
                    [GM_f,GM_g] = getGeometricMean(fxy_matrix,gxy_matrix,k1,k2);
                case 0
                    GM_f = 1;
                    GM_g = 1;
                    
            end
            
                    % Obtain normalised fxy and gxy matrices where
                    % binomials are included.
                    fxy_matrix_bi_n = fxy_matrix_bi./GM_f;
                    gxy_matrix_bi_n = gxy_matrix_bi./GM_g;
                    
                    % Obtain normalised fxy and gxy matrices.
                    fxy_matrix_n = fxy_matrix./GM_f;
                    gxy_matrix_n = gxy_matrix./GM_g;
            

            switch bool_thetas
                case 1
                    % Get the maximum and minimum entries of each
                    % coefficient of f, in the Sylvester matrix
                    [max_mtrx_f,min_mtrx_f] = GetMaxMin(fxy_matrix_n,n1,n2,k1,k2);

                    % Get the maximum and minimum entries of each
                    % coefficient of g, in the Sylvester matrix
                    [max_mtrx_g,min_mtrx_g] = GetMaxMin(gxy_matrix_n,n1,n2,k1,k2);

                    % Get an optimal value of alpha alone
                    %opt_alpha = OptimalAlpha(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g)
                    
                    % Get optimal values of theta1 and theta2 alone
                    [opt_theta_1, opt_theta_2] = OptimalTheta(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g);
                    
                    % Get optimal values of theta1 and theta2 and alpha
                    % [opt_alpha, opt_theta_1, opt_theta_2] = OptimalAlphaTheta(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g)
                    
                    % Add the optimal values of theta 1 and theta 2 to a
                    % matrix defined by k1,k2
                   
                case 0 % Exclude preprocessors
                    % set values of theta to 1.
                    opt_theta_1 = 1;
                    opt_theta_2 = 1;
            end
                    opt_theta_1_mtrx(k1+1,k2+1) = opt_theta_1;
                    opt_theta_2_mtrx(k1+1,k2+1) = opt_theta_2;
            
            Sk = BuildSubresultant(...
                    fxy_matrix_n, gxy_matrix_n,...
                    k1, k2,...
                    opt_theta_1_mtrx(k1+1,k2+1),...
                    opt_theta_2_mtrx(k1+1,k2+1));
            
            
            
            % add the minimum singular value S_{k1,k2} to the vector of
            % singular values for all k = k1+k2
            min_sing_val_vec(count) = min(svd(Sk));
            temp_vec_2(count,:) = [k1 k2];
            count = count + 1;
        end
    end
    
    % of all the minimal singular values for the given k1+k2 = k, stored in
    % the vector min_sing_val_vec, get the minimal singular value and its 
    % index.
    [sing_val_vec(k),index] = min(min_sing_val_vec);
    
    % get the values of k1 and k2 which gave the minimal value
    t_val_vec(k,:) = temp_vec_2(index,:);
    
end

opt_theta_1_mtrx;
opt_theta_2_mtrx;

% plot all the minimum singular values for k = 1,...,min(m,n)
plot(log10(sing_val_vec))

diff(sing_val_vec);

[~,maxindex] = max(diff(log10(sing_val_vec)));
degree_calc = maxindex;

[tval] = t_val_vec(degree_calc,:);
t1 = tval(1);
t2 = tval(2);


fprintf('The Calculated Degree of the GCD is given by \n %i \n \n',degree_calc)
fprintf('t1 = %i\n',t1)
fprintf('t2 = %i\n',t2)



% given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized

St = BuildSubresultant(fxy_matrix,gxy_matrix,t1,t2,opt_theta_1_mtrx(t1+1,t2+1),opt_theta_2_mtrx(t1+1,t2+1));

[rows,cols] = size(St);
% QR Decomposition of the Sylvester Matrix S_{k}
[Qk,Rk] = qr(St);
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
    ]  ;

n1-t1
n2-t2
m1-t1
m2-t2



% get coefficients of u and v
num_coeff_v = (n1-t1+1) * (n2-t2+1);
num_coeff_u = (m1-t1+1) * (m2-t2+1);
v_binom_calc = vecx(1:num_coeff_v);
u_binom_calc = vecx(num_coeff_v+1:end);

v_binom_calc = v_binom_calc./v_binom_calc(1)
u_binom_calc = u_binom_calc./u_binom_calc(1)

uxy_matrix./uxy_matrix(1,1)
vxy_matrix./vxy_matrix(1,1)



end




