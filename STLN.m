function [ fxy_output,gxy_output,X_output] = ...
    STLN( fxy_matrix,gxy_matrix,t1,t2,opt_col)
% Obtain the low rank approximation by method of SNTLN (Structured
% Non-Linear Total Least Norm) of the Sylvester matrix S(f,g).
%
%                             Inputs.
%
%
% fxy_matrix :    Coefficients of polynomial f, in standard bernstein basis.
%
% gxy_matrix :    Coefficients of polynomial g, in standard bernstein basis.
%
% t1 :   Degree of AGCD with respect to x.
%
% t2 :   Degree of AGCD with respect to y.
%
% lambda : 
%
% mu :
%
% opt_col : Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.
%
%
%                           Outputs.
%
%
% fx_output :- Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.
%
% gx_output :- Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.
%
%
% X_output :-
%
%

%                       Global Inputs

global MAX_ERROR_SNTLN
global MAX_ITERATIONS_SNTLN

if isempty(MAX_ERROR_SNTLN) || isempty(MAX_ITERATIONS_SNTLN)
    error('err')
end

%%

% Set the initial iterations number
ite = 1;

% Get degree of polynomials f.
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get degree of polynomial g.
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

% Get the number of coefficients in the polynomial f(x,y)
num_coeff_f = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
num_coeff_g = (n1+1) * (n2+1);

% Get the number of coefficients in both f(x,y) and g(x,y)
num_coeff = num_coeff_f + num_coeff_g;

% Get the number of coefficients in v(x,y)
num_coeff_v = (n1-t1+1) * (n2-t2+1);

% Get the number of coefficients in u(x,y)
num_coeff_u = (m1-t1+1) * (m2-t2+1);

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
num_coeff_x = num_coeff_u + num_coeff_v - 1;

% Create the identity matrix I, the matrix M formed from I by removing the
% column equivalent to the optimal column for removal from the Sylvester
% subresultant matrix, so that S_{t}(f,g)*M = A_{t}, where A_{t} is the
% Sylvester subresultant matrix with the column removed.

% Get the number of columns in C_{t}(f), the first partition of the Sylvester
% Matrix S_{t}(f,g)
num_cols_in_f_partition = (n1-t1+1) * (n2-t2+1);

% Get the number of columns in C_{t}(g), the second partition of the
% Sylvester matrix S_{t}(f,g)
num_cols_in_g_partition = (m1-t1+1) * (m2-t2+1);

% Get the total number of columns in the Sylvester matrix S_{t}(f,g)
num_cols = num_cols_in_f_partition + num_cols_in_g_partition;

% Create the identity matrix
I = eye(num_cols, num_cols);

% Create the matrix M, such that S(f,g)*M gives A_{t}, the Sylvester Matrix
% with the optimal column removed.
M = I;
M(:,opt_col) = [];

% Let e be the column removed from the identity matrix, such that
% S_{t}(f,g) * e gives the column c_{t}, where c_{t} is the optimal column
% removed from the Sylvester subresultant.
e = I(:,opt_col);



%% Form the Coefficient Matrix T = [C(f)|C(g)] such that DTQ * x = [col]

D = BuildD(t1,t2,m1,m2,n1,n2);

Q = BuildQ(n1-t1,n2-t2,m1-t1,m2-t2);

T1 = BuildT1(fxy_matrix,n1-t1,n2-t2);
T2 = BuildT1(gxy_matrix,m1-t1,m2-t2);
DTQ = D*[T1 T2]*Q;


%%
% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.

zk = zeros(num_coeff , 1);

% Get the coefficients of f(x,y) in matrix form
fww_vec = GetAsVector(fxy_matrix);
gww_vec = GetAsVector(gxy_matrix);

% Get the matrix A_{k}(f,g), which is the subresultant matrix S(f,g) with
% an opitmal column removed
At = DTQ;
ct = DTQ(:,opt_col);
At(:,opt_col) = [];

%%
% % Build the matrix of binomials corresponding to polynomial f(x,y)

% Multiply each row by nchoosek(m1,i)
% Multiply each column by nchoosek(m2,j)
Bi_m1 = GetBinomials(m1);
Bi_m2 = GetBinomials(m2);

binoms_f = diag(Bi_m1) * ones(m1+1,m2+1) * diag(Bi_m2);

% % Build the matrix of binomials corresponding to polynomial g(x,y)

% Multiply each row by nchoosek(n1,i)
% Multiply each column by nchoosek(n2,j)
Bi_n1 = GetBinomials(n1);
Bi_n2 = GetBinomials(n2);

binoms_g = diag(Bi_n1) * ones(n1+1,n2+1) * diag(Bi_n2);

% % Get the Binomials of f and g as vectors
binoms_f_vec = GetAsVector(binoms_f);
binoms_g_vec = GetAsVector(binoms_g);
G = diag([binoms_f_vec;binoms_g_vec]);
%%

P = BuildP(m1,m2,n1,n2,1,1,1,t1,t2,opt_col);
P = D*P*G;

%%
% Perform QR decomposition of Ak to obtain the solution x
x_ls = SolveAx_b(At,ct);

first_part = x_ls(1:(opt_col-1));
second_part = x_ls(opt_col:end);
x = [first_part ; 0 ; second_part];

% Build Matrix Y, where Y(v,u)*[f;g] = S(f,g)*[u;v]
Y = BuildY(m1,m2,n1,n2,t1,t2,opt_col,x_ls,1,1,1);
DYG = D*Y*G;

test1 = DYG * [fww_vec;gww_vec];
test2 = DTQ * x;
test1./test2;

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ct - (DTQ*M*x_ls);

%% Get the matrix p, which will store all the perturbations returned from LSE file
num_entries = num_coeff_f...
    + num_coeff_g ...
    + num_coeff_x ...
    + 3;


% Set the initial value of vector p to be zero
p = zeros(num_coeff_f...
    + num_coeff_g ...
    + num_coeff_x ...
    ,1);

%%
% Set the intial value of E to the identity matrix
E = eye(num_entries);

%%
% Create the matrix D(T+N)Q, initially N is empty so this is the same as T.

T1 = BuildT1(fxy_matrix,n1-t1,n2-t2);
T2 = BuildT1(gxy_matrix,m1-t1,m2-t2);

DTNQ = D*[T1 T2]*Q;

% Create the matrix C for input into iteration

H_z     = DYG-P;

H_x     = DTNQ*M;

C       = [H_z H_x];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    zk;...
    x_ls;...
    ];

yy              =   start_point;

% Set the termination criterion to a large value. It will be
% over written later.
condition(ite) = norm(res_vec)/norm(ct);

 % Edit 17/11/2015

xk = x_ls;

while condition(ite) >(MAX_ERROR_SNTLN) &&  ite < MAX_ITERATIONS_SNTLN
    % Use the QR decomposition to solve the LSE problem
    % min |y-p| subject to Cy=q
    
    y = LSE(E,p,C,res_vec);
    
    y_old = y;
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % obtain the small changes
    
    %% Break down y into its sections
    
    % get the coefficients corresponding to f and g
    delta_zk        = y(1:num_coeff_f + num_coeff_g ,1);
    % Remove the zk coefficients from the list of coefficients
    y(1:num_coeff_f + num_coeff_g) = [];
    
    % Get the coefficients corresponding to x
    delta_xk        = y(1:num_coeff_x,1);
    % Remove them from the list of coefficients
    y(1:num_coeff_x) = [];
    

    
    %% Update the variables
    
    % Update variables z_{k}, where z_{k} are perturbations in the
    % coefficients of f and g.
    zk = zk + delta_zk;
    
    % Update x_{k}, where x_{k} is the solution vector, containing
    % coefficients u and v.
    xk = xk + delta_xk;

    % Construct the Sylvester subresultant matrix DTQ.
    T1 = BuildT1(fxy_matrix,n1-t1,n2-t2);
    T2 = BuildT1(gxy_matrix,m1-t1,m2-t2);
    DTQ = D*[T1 T2]*Q;
    


    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ct = DTQ*e;
    

    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fx      = zk(1:num_coeff_f);
    z_gx      = zk(num_coeff_f + 1 :end);
    
    % Get the vectors z_fx and z_gx as matrices, which match the shape of
    % f(x) and g(x).
    z_fx_mat = GetAsMatrix(z_fx,m1,m2);
    z_gx_mat = GetAsMatrix(z_gx,n1,n2);
    
    
    % Build the coefficient Matrix N = [T(z1) T(z2)], of structured perturbations, with
    % same structure as DTQ.
    N1 = BuildT1(z_fx_mat,n1-t1,n2-t2);
    N2 = BuildT1(z_gx_mat,m1-t1,m2-t2);
    DNQ = D*[N1 N2]*Q;

    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{k} the column of S_{k} moved to the RHS
    h = DNQ*e;
       
    % Build the matrix (T+N)
    TN1 = BuildT1(fxy_matrix + z_fx_mat,n1-t1,n2-t2);
    TN2 = BuildT1(gxy_matrix + z_gx_mat,m1-t1,m2-t2);
    DTNQ = D*[TN1 TN2]*Q;

    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    Y = BuildY(m1,m2,n1,n2,t1,t2,opt_col,xk,1,1,1);
    
    DYG = D*Y*G;
    
    test1 = DYG * [fww_vec;gww_vec];
    first_part = xk(1:(opt_col-1));
    second_part = xk(opt_col:end);
    x = [first_part ; 0 ; second_part];
    test2 = DTQ * x;
    
    
    % Calculate the matrix DP where P is the matrix such that c = P[f;g]
    P = BuildP(m1,m2,n1,n2,1,1,1,t1,t2,opt_col);
    P = D*P*G;
    
    % Get residual as a vector
    rk = (ct+h) - DTNQ*M*xk ;
    
    % Create the matrix C. This is made up of five submatrices, HZ, Hx,
    % H_alpha and H_theta1 and H_theta2.
    
    Hz          = DYG-P;
    
    Hx          = DTNQ*M;
    
   
    
    C = [Hz,Hx];  % the matrix C
    
    % Calculate the new right hand vector
    ek = ct + h;
    
    % Update gnew - used in lse problem
    res_vec = rk;
    
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(rk) / norm(ek);
    
    % Update fnew - used in LSE Problem.
    p = -(yy-start_point);
    
    
end

% Plot Graphs
global PLOT_GRAPHS
switch PLOT_GRAPHS
    case 'y'
        figure('name','Residuals in SNTLN')
        hold on
        title('Residuals in SNTLN')
        xlabel('Iterations')
        ylabel('log_{10} Residuals')
        plot((1:1:ite),log10(condition),'-s')
        hold off
        
       
    case 'n'
    otherwise
        error('err')
end

%%
if ite == MAX_ITERATIONS_SNTLN
    fprintf('SNTLN Failed to converge, default to input values\n')
    fxy_output = fxy_matrix;
    gxy_output = gxy_matrix;
    X_output = x_ls;
    return;
end

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.

% get the vector zk

zPert_f_vec = zk(1:num_coeff_f);
zPert_f_mat = GetAsMatrix(zPert_f_vec,m1,m2);

zPert_g_vec = zk(num_coeff_f+1:end);
zPert_g_mat = GetAsMatrix(zPert_g_vec,n1,n2);

fxy_output = fxy_matrix + zPert_f_mat;

gxy_output = gxy_matrix + zPert_g_mat;

X_output  = xk;


% Print the number of iterations
fprintf('--------------------------------------------------------------------------- \n')
fprintf('Iterations over Sylvester Matrix : %i \n', ite);
fprintf('--------------------------------------------------------------------------- \n')
end

















