function [] = o_gcd_new_new(ex_num)

global bool_geomean
bool_geomean = 0;

global bool_method

% bool_method = 0 - Include all columns in the subresultant matrices.
bool_method = 1;

% Get the example file
[f_x_roots,f_y_roots,...
    g_x_roots,g_y_roots,...
    u_x_roots, u_y_roots,...
    v_x_roots, v_y_roots,...
    d_x_roots, d_y_roots,...
    m1,m2,n1,n2,t1,t2] = examples(ex_num);

% Get the total degrees
M = m1+m2;
N = n1+n2;
T_exact = t1+t2;

% Get total degrees of polynomial u and v
% Get total degree of polynomial u
M_T = m1-t1 + m2-t2;
% Get total degree of polynomial v
N_T = n1-t1 + n2-t2;

% 2.1 get coefficients of polynomial f
% % 2.1.1 get coefficients of polynomial f in terms of x
% % 2.1.2 get coefficients of polynomial f in terms of y
% % 2.1.3 Build matrix of coefficients of f(xy)

f_x_poly = BuildPoly(f_x_roots);
f_y_poly = BuildPoly(f_y_roots);
fxy_matrix = f_y_poly' * f_x_poly
fxy_vec = [];
for i = -m2 : 1 : m1
    temp_diag = flipud(diag(flipud(fxy_matrix),i));
    fxy_vec = [fxy_vec ; temp_diag];
end

% 2.2 get coefficients of polynomial g
% % 2.2.1 get coefficients of polynomial g in terms of x
% % 2.2.2 get coefficients of polynomial g in terms of y
% % 2.2.3 get coefficients of polynomial g in terms of x and y
g_x_poly = BuildPoly(g_x_roots);
g_y_poly = BuildPoly(g_y_roots);
gxy_matrix = g_y_poly' * g_x_poly;

gxy_vec = [];
for i = -n2 : 1 : n1
    temp_diag = flipud(diag(flipud(gxy_matrix),i));
    gxy_vec = [gxy_vec ; temp_diag ];
end

% 2.3 get coefficients of polynomial u
% % 2.3.1 get coefficients of polynomial u in terms of x
% % 2.3.2 get coefficients of polynomial u in terms of y
% % 2.3.3 get matrix of coefficients of u in terms of x and y

u_x_poly = BuildPoly(u_x_roots);
u_y_poly = BuildPoly(u_y_roots);
uxy_matrix = u_y_poly' * u_x_poly;
uxy_vec = [];
try
    for i = -(m2-t2) : 1 : (m1-t1)
        flipud(diag(flipud(uxy_matrix),i));
        uxy_vec = [uxy_vec ; flipud(diag(flipud(uxy_matrix),i))];
    end
catch
    uxy_vec = uxy_matrix;
end

% 2.4 get coefficients of polynomial v
% % 2.4.1 get coefficients of polynomial v in terms of x
% % 2.4.2 get coefficients of polynomial v in terms of y
% % 2.4.3 get coefficients of polynomial v in terms of x and y
v_x_poly = BuildPoly(v_x_roots);
v_y_poly = BuildPoly(v_y_roots);
vxy_matrix = v_y_poly' * v_x_poly;
vxy_vec = [];

try
    for i = -(n2-t2) : 1 : (n1-t1)
        i;
        flipud(diag(flipud(vxy_matrix),i));
        vxy_vec = [vxy_vec ; flipud(diag(flipud(vxy_matrix),i))];
    end
catch
    vxy_vec = vxy_matrix;
end

% % Initalise an empty vector Y to store values from QR Decomposition

Y = [];





% For each total degree k = 1,...,min(Tm,Tn)
for k = 1:1:min(M,N)
    
    % Build Sylvester Matrix
    Cauchy_f = BuildCauchy(fxy_matrix,m1,m2,n1,n2,k);
    Cauchy_g = BuildCauchy(gxy_matrix,n1,n2,m1,m2,k);
    
    switch bool_geomean
        case 1
            % Normalise coefficients by geometric mean
            Cauchy_f = Cauchy_f./geomean(abs(Cauchy_f(Cauchy_f~=0)));
            Cauchy_g = Cauchy_g./geomean(abs(Cauchy_g(Cauchy_g~=0)));
    end
    
    rows = max([size(Cauchy_f,1) size(Cauchy_g,1)]);
    colsf = size(Cauchy_f,2);
    colsg = size(Cauchy_g,2);
    
    padded_f = zeros(rows,colsf);
    padded_f(1:size(Cauchy_f,1),:) = Cauchy_f;
    
    padded_g = zeros(rows,colsg);
    padded_g(1:size(Cauchy_g,1),:) = Cauchy_g;
    
    S_k = [padded_f, padded_g];
    
    
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
plot(log10(min_sing_value),'black-s');

hold off

% % Get degree of GCD
[~,degree_calc_maxmin_diag] = max(abs(diff(log10(max_over_min_diag)))) ;
degree_calc_maxmin_diag = degree_calc_maxmin_diag ;
fprintf('Max/Min Ratio method to obtain degree: \n t = %i \n',degree_calc_maxmin_diag);

% Set the degree t
T = degree_calc_maxmin_diag;

fprintf('Total Degree Exact:\n')
T_exact
fprintf('Total Degree Calc:\n')
degree_calc_maxmin_diag

fprintf('Set t to be exact \n')


T = T_exact;


% Given the calculated degree, build the subresultant and obtain the
% polynomials u and v from it.
Cauchy_f = BuildCauchy(fxy_matrix,m1,m2,n1,n2,T);
Cauchy_g = BuildCauchy(gxy_matrix,n1,n2,m1,m2,T);

rows = max([size(Cauchy_f,1) size(Cauchy_g,1)])
colsf = size(Cauchy_f,2);
colsg = size(Cauchy_g,2);

padded_f = zeros(rows,colsf);
padded_f(1:size(Cauchy_f,1),:) = Cauchy_f;

padded_g = zeros(rows,colsg);
padded_g(1:size(Cauchy_g,1),:) = Cauchy_g;

S_t = [padded_f, padded_g]

% Get optimal column for removal from S_{k}
% 2.3 Remove each column of the subresultant in turn
res_arr = [];
for q = 1:1:size(S_t,2)
    
    % 2.3.1 Obtain column to be removed
    ck = S_t(:,q) ;
    
    % 2.3.2 Copy sk to Ak and remove the optimal column
    A_t = S_t;
    A_t(:,q) = [];
    
    % 2.3.3 Obtain x, the approximate solution vector
    x = pinv(A_t)*ck;
    
    % 2.3.4 Normalise the residual
    res_arr(q)=norm(ck-(A_t*x));
end

% 2.4 Get column for which residual was minimal
[~,opt_col] = min(log10(res_arr));


% 2.5 Calculate the quotient polynomials u and v
% Obtain column to be removed
ck = S_t(:,opt_col) ;

% 2.3.2 Copy sk to Ak and remove the optimal column
A_t = S_t;
A_t(:,opt_col) = [];

% Obtain least squares solution
x_ls = pinv(A_t) * ck;

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ]  

% Obtain values for quotient polynomials u and v. still expressed in the
% scaled bernstein basis, including theta.




switch bool_method
    case 0 % Use all coefficients
        vx = vecx(1:(N-T+1)*(N-T+1));
        ux = -vecx((N-T+1)*(N-T+1)+1:end);
    case 1 % Use trimmed coefficients
        vx = vecx(1:nchoosek(N-T+2,N-T));
        ux = -vecx(nchoosek(N-T+2,N-T)+1:end);
end

fprintf('Coefficients of v(x,y) Calculated \n')
vx./vx(1)
fprintf('Coefficients of v(x,y) Exact \n')
vxy_vec./vxy_vec(1)
vxy_matrix;
fprintf('Coefficients of u(x,y) Calculated \n')
ux./ux(1)

fprintf('Coefficients of u(x,y) Exact \n')
uxy_vec./uxy_vec(1)
uxy_matrix;
% Obtain coprime polynomials u and v




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

function [Cauchy_f] = BuildCauchy(fxy_matrix,m1,m2,n1,n2,t)


% Make use of global variable
 global bool_method

 
 % where t is the index

% Get total degree of polynomial g
T_degree_g = n1+n2;

% Get total degree of GCD(f,g)
T_degree_d = t;

% Get total degree of polynomial v
T_degree_v = T_degree_g - T_degree_d;

% Build the matrix such that the coefficients of f when multiplied by the
% basis of v, fill the entries.
 fxy_matrix_zeros = zeros((m2+1)+T_degree_v,(m1+1)+T_degree_v);
%fxy_matrix_zeros = zeros((m1)+T_degree_v,(m1+1)+T_degree_v);

fxy_matrix_zeros(1:size(fxy_matrix,1),1:size(fxy_matrix,2)) = fxy_matrix;

% initialise a cell array to store all the matrices of f*basis(v)
mycellarray = cell(T_degree_v+1,T_degree_v+1);

% initialise a vector array to store all the vectors of f*basis(v) in
% vector format.
myvecarray = cell(T_degree_v+1,T_degree_v+1);

Cauchy_f = [];

% for each column j in vxy matrix
for j = 0:1:T_degree_v
    % for each row i in vxy matrix
    for i = 0:1:T_degree_v
        
        temp_matrix = circshift(fxy_matrix_zeros,[i,j]);
        mycellarray{i+1,j+1} = temp_matrix;
        
        % Get vector of antidiagonals of temp matrix
        temp_vec = [];
        for k = -(T_degree_v+1)-1-1-10 : 1 : T_degree_v+1+1+1+10;
            temp_diag = flipud(diag(flipud(temp_matrix),k));
            temp_vec = [temp_vec ; temp_diag];
        end
        temp_vec;
        myvecarray{i+1,j+1} = temp_vec;
        
    end
end


Cauchy_f;

% of the 36 values only the first 21 are included, that is the values in
% the upper triangle
included_cols = nchoosek(T_degree_v+2,T_degree_v);


for k = 2:1:((T_degree_v+1)  + (T_degree_v+1))
    for i = 1:1:k
        j = k-i;
        if i > T_degree_v+1
        elseif j > T_degree_v+1
        elseif j < 1
            break;
        else
            temp_vec = myvecarray{i,j};
            Cauchy_f = [Cauchy_f temp_vec];
        end
    end
    
end


switch bool_method
    case 0 % include all columns
        
    case 1 % only include the nonzero columns
        Cauchy_f = Cauchy_f(:,1:included_cols);
end
end


function [t]=B_conv(rooot,mult)
%% This function convolves the vector [-r 1-r] with itself m times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% r : root
% m : multiplicity of root
% Outputs:
% t : vector which stores the result from this convolution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Note that (y-r) = -r(1-y) + (1-r)y and thus the polynomial y-r in the
% power basis is represented as the polynomial -r(1-y) + (1-r)y in the
% scaled Bernstein basis.


if mult==1
    t=[-rooot,1-rooot];
else
    
    q=[-rooot,1-rooot];
    t=[-rooot,1-rooot];
    for k=2:1:mult
        t=conv(t,q);
    end
end
end
