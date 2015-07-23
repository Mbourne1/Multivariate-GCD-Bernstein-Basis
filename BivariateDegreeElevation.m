function [] = BivariateDegreeElevation(ex_num)

% % Variables Used in this file

% f_x_roots - Roots of f in the x variable

% f_y_roots - Roots of f in the y variable

% g_x_roots - Roots of g in the x variable

% g_y_roots - Roots of g in the y variable

% d_x_roots - Roots of d in the x variable

% d_y_roots - Roots of d in the y variable

% fxy_matrix    - Matrix of coefficients of fxy where the rows represent
% the basis elements 1,x,...x^m1, and the column 1,y,...,y^m2

% fxy_matrix_bi - Matrix fxy_matrix with the corresponding binomial
% coefficients included.

% fxy_binoms    - Matrix of binomial coefficients corresponding to the
% coefficiens of fxy_matrix.

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
fxy_matrix_bi =  f_x_poly' * f_y_poly;

% for each row of fxy_matrix_bi. degree elevate

% get matrix of polynomial g coefficients including the binomial
% coefficients.
gxy_matrix_bi =  g_x_poly' * g_y_poly;

% get matrix of polynomial coefficients including the binomial coefficients
uxy_matrix_bi = u_x_poly' * u_y_poly;
vxy_matrix_bi = v_x_poly' * v_y_poly;

% we must degree elevate each row of the fxy_matrix so that it forms an
% upper triangular matrix


% build the leading diagonal matrix of basis entries in terms of x

[r,c] = size(fxy_matrix_bi)

getInPowerBasis(fxy_matrix_bi)
gxy_matrix_bi
uxy_matrix_bi
vxy_matrix_bi


fprintf('Degree of input polynomial f: m1 = %i, m2 = %i, m = %i\n' ,m1,m2,m)
fprintf('Degree of input polynomial g: n1 = %i, n2 = %i, n = %i\n' ,n1,n2,n)
fprintf('Degree of polynomial d: d1 = %i, d2 = %i, d = %i \n' ,d1,d2,d)

% % Build a matrix which is equal only to the binomial coefficients of
% fxy_matrix

% . Initialise the Empty matrix
fxy_binoms = zeros(m1+1,m2+1);
% . For each row of the matrix fxy - for each basis element 1,x,...x^m1
for i = 0:1:m1
    % . For each column of the matrix fxy - for each basis element
    % 1,y,...,y^m2
    for j = 0:1:m2
        fxy_binoms(i+1,j+1) = nchoosek(m1,i) * nchoosek(m2,j);
    end
end



% % Build a matrix which is equal only to the binomial coefficiens of
% gxy_matrix

% . Initialise the empty matrix
gxy_binoms = zeros((n1+1),(n2+1));
% . For each row of the matrix gxy - for each basis element 1,x,...,x^m1
for i = 0:1:n1
    % . For each column of the matrix gxy - for each basis element
    % 1,y,...,y^m2
    for j = 0:1:n2
        gxy_binoms(i+1,j+1) = nchoosek(n1,i) * nchoosek(n2,j);
    end
end

fprintf('The binomial matrices')
fxy_binoms
gxy_binoms

% Strip fxy_matrix_bi of binomial coefficients to obtain fxy_matrix
fxy_matrix = fxy_matrix_bi./fxy_binoms;

% Strip gxy_matrix_bi of binomial coefficients to obtain fxy_matrix
gxy_matrix = gxy_matrix_bi./gxy_binoms;

fprintf('f and g with binomials removed')
fxy_matrix
gxy_matrix

% Get number of additional entries in x and y for degree elevation

        deg_elv_fxy_matrix_bi = fxy_matrix_bi;
        deg_elv_gxy_matrix_bi = gxy_matrix_bi;



Y = [];
%
% Let K be defined to be t
%
for k = 1:1:min(m,n)
    
    % Build Matrices to cast fxy matrix and gxy matrix into. such that they
    % are surrounded by zeros
    
    cast_fxy_bi = zeros(m+n-k+1,m+n-k+1);
    cast_gxy_bi = zeros(m+n-k+1,m+n-k+1);
    
    cast_fxy_bi(1:size(deg_elv_fxy_matrix_bi,1),1:size(deg_elv_fxy_matrix_bi,2)) =...
        deg_elv_fxy_matrix_bi;
    
    cast_gxy_bi(1:size(deg_elv_gxy_matrix_bi,1),1:size(deg_elv_gxy_matrix_bi,2)) =...
        deg_elv_gxy_matrix_bi;
    
    
    % Build the cauchy matrix for polynomial f
    
    Cauchy_f = BuildCauchy(n,m,k,cast_fxy_bi)
    Cauchy_g = BuildCauchy(m,n,k,cast_gxy_bi);
    
    % Build Matrix D^{-1} the diagonal matrix
    Bool_D = 1;
    switch Bool_D
        case 0
            D = eye(nchoosek(m+n-k+2,2),nchoosek(m+n-k+2,2));
        case 1
            X = [];
            temp_vec = [];
            for k2 = 0:1:m+n-k
                for i = k2 :-1:0
                    j = k2-i;
                    X = [X ;[i,j]];
                    temp_vec = [temp_vec ; nchoosek(m+n-k,i) * nchoosek(m+n-k,j)];
                end
            end
            D = diag(1./temp_vec);
    end
    
    % Build the Sylvester matrix
    
    if k == 3
    size(Cauchy_f)
        Cauchy_f
    
    end
    
    
    S_k = D*[Cauchy_f Cauchy_g];
    
    figure(999)
    hold on
    plot(log10(svd(S_k)))
    hold off
    
    % Obtain Singular Value decomposition
    min_sing_value(k) = min(svd(S_k));
    
    % For each subresultant perform QR decomposition
    [Q,R] = qr(S_k);
    
    % Get absolute values of R
    R = abs(R);
    
    % Get number of rows in R1
    [R1_rows,~] = size(diag(R));
    
    % Obtain R1
    R1 = R(1:R1_rows,1:R1_rows);
    
    % Normalise R1
    R1 = R1./norm(R1);
    R1 = R1./R1(1);
    % Obtain sum of the rows of R1
    R1_sum_rows = sum(R1,2);
    
    % Get the diagonal elements of rows of R1
    R1_diag = abs(diag(R));
    
    % Normalise all elements of the matrix R
    Rs = sqrt(sum(abs(R1_sum_rows).^2,2));
    
    % Scatter plot for singular values
    
    %Scatter Plot Data
    ks = k.*ones(size(Rs));
    ns = 1:1:size(Rs,1);
    X = [ks Rs ns'];
    Y = [Y; X];
    
    max_over_min_diag(k) = max(R1_diag)./min(R1_diag);
    
    % for each column in the sylvester subresultant matrix, remove it to RHS
    for q = 1:1:size(S_k,2)
        
        % Obtain column to be removed
        ck = S_k(:,q) ;
        
        % Copy sk to Ak and remove the optimal column
        Ak = S_k;
        Ak(:,q) = [];
        
        % Obtain x, the approximate solution vector
        x = pinv(Ak)*ck;
        
        % Normalise the residual
        res_arr(q)=norm(ck-(Ak*x));
    end
    
    % get minimal residual from the given subresultant
    residual_arr(k) = min(res_arr);
    
end


figure(2)
hold on
scatter(Y(:,1),log10(Y(:,2)))
hold off

figure(3)
hold on
title('Residuals obtained by removing optimal column from each subresultant S_{k}')
xlabel('k = index of Subresultant S_{k}')
ylabel('log_{10} residual')
plot(log10(residual_arr),'-s')
hold off

figure(4)
hold on
title('max:min diag elements')
plot(log10(max_over_min_diag),'red-s');
hleg = legend('max/min diag element');
title('diag elements');
hold off

figure(5)
hold on
title('minimum singular values of S_{k}')
plot((min_sing_value),'black-s');
hold off

% % Get degree of GCD
[~,degree_calc_maxmin_diag] = max(abs(diff(log10(max_over_min_diag)))) ;
degree_calc_maxmin_diag = degree_calc_maxmin_diag ;
fprintf('Max/Min Ratio method to obtain degree: \n t = %i \n',degree_calc_maxmin_diag);

% Set the degree t
T = degree_calc_maxmin_diag;

fprintf('Total Degree Exact: %i \n', d)
d;
fprintf('Total Degree Calc: %i \n',degree_calc_maxmin_diag )
degree_calc_maxmin_diag;

t = degree_calc_maxmin_diag ;

% Build the Tth subresultant matrix where t is the calculated degree of the
% GCD.

cast_fxy_bi = zeros(m+n-t+1,m+n-t+1);
cast_gxy_bi = zeros(m+n-t+1,m+n-t+1);

cast_fxy_bi(1:size(deg_elv_fxy_matrix_bi,1),1:size(deg_elv_fxy_matrix_bi,2)) =...
    deg_elv_fxy_matrix_bi;

cast_gxy_bi(1:size(deg_elv_gxy_matrix_bi,1),1:size(deg_elv_gxy_matrix_bi,2)) =...
    deg_elv_gxy_matrix_bi;

% Build the cauchy matrix for polynomial f

Cauchy_f = BuildCauchy(n,m,t,cast_fxy_bi);
Cauchy_g = BuildCauchy(m,n,t,cast_gxy_bi);

% Build Matrix D^{-1} the diagonal matrix
Bool_D = 1;
switch Bool_D
    case 0
        D = eye(nchoosek(m+n-t+2,2),nchoosek(m+n-t+2,2));
    case 1
        X = [];
        temp_vec = [];
        for k2 = 0:1:m+n-t
            for i = k2 :-1:0
                j = k2-i;
                X = [X ;[i,j]];
                temp_vec = [temp_vec ; nchoosek(m+n-t,i) * nchoosek(m+n-t,j)];
            end
        end
        D = diag(1./temp_vec);
end

% Build the t-th subresultant matrix where t is the calculated degree
S_t = D*[Cauchy_f Cauchy_g];

% for each column in the sylvester subresultant matrix, remove it to RHS
for q = 1:1:size(S_t,2)
    
    % Obtain column to be removed
    ck = S_t(:,q) ;
    
    % Copy sk to Ak and remove the optimal column
    At = S_t;
    At(:,q) = [];
    
    % Obtain x, the approximate solution vector
    x = pinv(At)*ck;
    
    % Normalise the residual
    res_arr(q)=norm(ck-(At*x));
end

% get minimal residual from the given subresultant
residual_arr(k) = min(res_arr);

% Get column of minimal residual
[~,opt_col] = min(residual_arr);

% Obtain column to be removed
ck = S_t(:,q) ;

% Copy sk to Ak and remove the optimal column
A_k = S_t;
A_k(:,q) = [];

% Obtain x, the approximate solution vector
x = pinv(A_k)*ck;

% split x
vecx =[
    x(1:(opt_col)-1);
    -1;
    x(opt_col:end);
    ]  ;

% polynomial v has \binom{n-k+2}{2} coefficients
% polynomial u has \binom{m-k+2}{2} coefficients

v_bi = vecx(1:nchoosek(n-t+2,2),:);
u_bi = -vecx(nchoosek(n-t+2,2)+1:end,:);

fprintf('Calculated u and v in vector form\n')
fprintf('Coefficients of u\n')
fprintf('%30.15f \n', [u_bi./u_bi(1)])
fprintf('Coefficients of v\n')
fprintf('%30.15f \n', [v_bi./v_bi(1)])




% Normalise the exact uxy_matrix and vxy_matrix
fprintf('Normalised uxy_matrix and vxy_matrix (Exact)')
uxy_matrix_bi./(uxy_matrix_bi(1))
vxy_matrix_bi./(vxy_matrix_bi(1))

getInPowerBasis(uxy_matrix_bi)

getInPowerBasis(vxy_matrix_bi)

end


function [f_bi] = BuildPoly(A)

% Calculate the number of distinct roots of the polynomial.
r = size(A,1);

% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.

f_bi = 1;
for k = 1:1:r
    w = B_conv(A(k,1),A(k,2));
    f_bi = conv(f_bi,w) ;
end


end


function sum = summing_function(fxy_matrix,m,n,r,s,i,j)

sum = 0;
for k = max(0,i-r) : 1 : min(m,i)
    for l = max(0,j-s) : 1 : min(n,j)
        sum = sum +...
            (...
            fxy_matrix(k+1,l+1) *...
            nchoosek(m,k) *...
            nchoosek(r,i-k) *...
            nchoosek(n,l) *...
            nchoosek(s,j-l) /...
            ...
            (nchoosek(m+r,i) *...
            nchoosek(n+s,j))...
            );
    end
end

end

function Cauchy_f = BuildCauchy(m,n,k,cast_fxy_bi)

Cauchy_f = [];

for p = 0:1:n-k
    for i = p:-1:0
        j = p-i;
        temp_matrix = circshift(cast_fxy_bi,[i,j]);
        temp_vec = [];
        %for z = -(m+n-k+1) : 1 : (m+n-k+1)
        for z = -(m+n-k+1) : 1 : (0)
            temp_diag = diag(flipud(temp_matrix),z);
            temp_vec = [temp_vec ; temp_diag];
        end
        Cauchy_f = [Cauchy_f temp_vec];
    end
end
end


function deg_elv_fxy_matrix = DegreeElevate(fxy_matrix,m1,m2,r_f,s_f)

m = m1+m2;

% Degree elevate fxy_matrix
deg_elv_fxy_matrix = zeros(m+1,m+1);
% for each column j
for j = 0:1:m2+s_f
    % for each row i
    for i = 0:1:m1+r_f
        deg_elv_fxy_matrix(i+1,j+1) = summing_function(fxy_matrix,m1,m2,r_f,s_f,i,j);
    end
end

end



