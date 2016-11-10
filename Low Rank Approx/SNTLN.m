function [ fxy_output,gxy_output,alpha_output,theta1_output,theta2_output,X_output] = ...
    SNTLN( fxy_matrix,gxy_matrix,...
    i_alpha, i_theta1, i_theta2,...
    k1,k2,...
    idx_col)
% Obtain the low rank approximation by method of SNTLN (Structured
% Non-Linear Total Least Norm) of the Sylvester matrix S(f,g).
%
% % Inputs.
%
% fx_matrix : Coefficients of polynomial f, in standard bernstein basis.
%
% gx_matrix : Coefficients of polynomial g, in standard bernstein basis.
%
% i_alpha :   Initial value of alpha
%
% i_theta :   Initial value of theta
%
% k1 : Degree of polynomial d(x,y) with respect to x
%
% k2 : Degree of polynomial d(x,y) with respect to y
%
% idx_col : Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.
%
% % Outputs.
%
% fx_output :- Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.
%
% gx_output :- Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.
%
% alpha_output :-
%
% theta_output :-
%
% X_output :-
%
%

% Global Inputs


global SETTINGS

if isempty(SETTINGS.MAX_ERROR_SNTLN) || isempty(SETTINGS.MAX_ITERATIONS_SNTLN)
    error('err')
end

%%

% Set the initial iterations number
ite = 1;

% Set initial values of alpha and theta
th1(ite) = i_theta1;
th2(ite) = i_theta2;
alpha(ite) = i_alpha;

% Get degree of polynomials f.
[m1,m2] = GetDegree(fxy_matrix);

% Get degree of polynomial g.
[n1,n2] = GetDegree(gxy_matrix);

% Get the number of coefficients in the polynomial f(x,y)
nCoeff_fxy = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
nCoeff_gxy = (n1+1) * (n2+1);

% Get the number of coefficients in both f(x,y) and g(x,y)
nCoeff_fg = nCoeff_fxy + nCoeff_gxy;

% Get the number of coefficients in v(x,y)
nCoeff_vxy = (n1-k1+1) * (n2-k2+1);

% Get the number of coefficients in u(x,y)
nCoeff_uxy = (m1-k1+1) * (m2-k2+1);

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
nCoeff_x = nCoeff_uxy + nCoeff_vxy - 1;

% Get the number of columns in C_{n-k}(f), the first partition of the Sylvester
% Matrix S_{t}(f,g).
nCols_Tf = nCoeff_vxy;

% Get the number of columns in C_{t}(g), the second partition of the
% Sylvester matrix S_{t}(f,g)
nCols_Tg = nCoeff_uxy;

% Get the total number of columns in the Sylvester matrix S_{t}(f,g)
nCols_Sk = nCols_Tf + nCols_Tg;

% Create the identity matrix
I = eye(nCols_Sk, nCols_Sk);

% Create the matrix M, such that S(f,g)*M gives A_{t}, the Sylvester Matrix
% with the optimal column removed.
M = I;
M(:,idx_col) = [];

% Let e be the column removed from the identity matrix, such that
% S_{t}(f,g) * e gives the column c_{t}, where c_{t} is the optimal column
% removed from the Sylvester subresultant.
e = I(:,idx_col);

%% Preprocessing

% Normalise polynomial f(x,y) by geometric means

% Obtain polynomials in Modified Bernstein Basis, using initial values of
% alpha and theta.

% Multiply the rows of fxy_matrix by theta1, and multiply the cols of
% fxy_matrix by \theta_{2}
fww_matrix = GetWithThetas(fxy_matrix_n,th1,th2);

% Multiply the rows of gxy_matrix by theta1, and multiply the cols of
% gxy_matrix by theta2
gww_matrix = GetWithThetas(gxy_matrix_n,th1,th2);

%% Form the Coefficient Matrix T = [C(f)|C(g)] such that DTQ * x = [col]

D = BuildD(k1,k2,m1,m2,n1,n2);
DTQ_fg = BuildDTQ(fww_matrix,alpha(ite).*gww_matrix,k1,k2);

% Get partial derivative of f(\omega_{1},\omega_{2}) with respect to \alpha
fww_wrt_alpha            = zeros(m1+1,m2+1);

% Get partial derivative of g(\omega_{1},\omega_{2}) with respect to \alpha
alpha_gww_wrt_alpha      = gxy_matrix;

% Calculate the partial derivatives of f(w,w) with respect to theta_1
fww_wrt_th1 = Differentiate_wrt_theta1(fww_matrix,th1(ite));

% Calculate the partial derivative of g(w,w) with respect to theta_1
gww_wrt_th1 = Differentiate_wrt_theta1(gww_matrix,th1(ite));

% Calculate the partial derivative of f(w,w) with respect to theta_2
fww_wrt_th2 = Differentiate_wrt_theta2(fww_matrix,th2(ite));

% Calculate the partial deriviates of g(w,w) with respect to theta_2
gww_wrt_th2 = Differentiate_wrt_theta2(gww_matrix,th2(ite));

% Build the derivative of T(f,g) with respect to alpha
DTQ_alpha = BuildDTQ(fww_wrt_alpha,alpha_gww_wrt_alpha,k1,k2);

% Calculate the derivative of T(f,g) with respect to theta_{1}
DTQ_wrt_th1 = BuildDTQ(fww_wrt_th1, alpha(ite).*gww_wrt_th1,k1,k2);

% Calcualte the derivative of T(f,g) with respect to theta_2
DTQ_wrt_th2 = BuildDTQ(fww_wrt_th2, alpha(ite).*gww_wrt_th2,k1,k2);


%%
% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.

nRows_Sylv_mat = (m1 + n1 - k1 + 1) * (m2 + n2 - k2 + 1);
nCols_T1 = (n1 - k1 + 1) * (n2 - k2 + 1);
nCols_T2 = (m1 - k1 + 1) * (m2 - k2 + 1);

zk = zeros(nCoeff_fg , 1);

%%
% Initilaise the derivative of N wrt alpha.
DNQ_wrt_alpha   = zeros(nRows_Sylv_mat,nCols_T1 + nCols_T2);

% Initilaise the derivative of N wrt theta_1.
DNQ_wrt_theta1   = zeros(nRows_Sylv_mat,nCols_T1 + nCols_T2);

% Initialise the derivative of N wrt theta 2
DNQ_wrt_theta2   = zeros(nRows_Sylv_mat,nCols_T1 + nCols_T2);

%%
% Initialise the derivative of h
% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
ht_wrt_alpha     = DNQ_wrt_alpha*e;
ht_wrt_theta1    = DNQ_wrt_theta1*e;
ht_wrt_theta2    = DNQ_wrt_theta2*e;

% Get the matrix A_{k}(f,g), which is the subresultant matrix S(f,g) with
% an opitmal column removed
Ak_fg = DTQ_fg;
ck = DTQ_fg(:,idx_col);
Ak_fg(:,idx_col) = [];

%%
% % Build the matrix of binomials corresponding to polynomial f(x,y)

binomials_f = GetWithBinomials(ones(m1+1,m2+1));
vBinomials_f = GetAsVector(binomials_f);

% % Build the matrix of binomials corresponding to polynomial g(x,y)
binomials_g = GetWithBinomials(ones(n1+1,n2+1));
vBinomials_g = GetAsVector(binomials_g);

% % Get the Binomials of f and g as vectors
G = diag([vBinomials_f;vBinomials_g]);

% %
P = BuildP(m1,m2,n1,n2,th1(ite),th2(ite),alpha(ite),k1,k2,idx_col);
DPG = D*P*G;

%%
% Calculate the derivatives wrt alpha and theta of the removed column.
ck_wrt_alpha        = DTQ_alpha*e;
ck_wrt_theta1       = DTQ_wrt_th1*e;
ck_wrt_theta2       = DTQ_wrt_th2*e;


%%
% Perform QR decomposition of Ak to obtain the solution x
x_ls = SolveAx_b(Ak_fg,ck);


% Build Matrix Y, where Y(v,u)*[f;g] = S(f,g)*[u;v]
Y = BuildY(m1,m2,n1,n2,k1,k2,idx_col,x_ls,alpha(ite),th1(ite),th2(ite));
DYG = D*Y*G;

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ck - (DTQ_fg*M*x_ls);

%% Get the matrix p, which will store all the perturbations returned from LSE file
nEntries = nCoeff_fxy...
    + nCoeff_gxy ...
    + nCoeff_x ...
    + 3;


% Set the initial value of vector p to be zero
p = zeros(nCoeff_fxy...
    + nCoeff_gxy ...
    + nCoeff_x ...
    + 3,1);

%%
% Set the intial value of E to the identity matrix
E = eye(nEntries);

% Create the matrix D(T+N)Q, initially N is empty so this is the same as T.
DTNQ = BuildDTQ(fww_matrix, alpha(ite).*gww_matrix,k1,k2);

% Create The matrix (T+N) with respect to alpha
DTNQ_wrt_alpha = BuildDTQ(fww_wrt_alpha, alpha_gww_wrt_alpha,k1,k2);

% Create The matrix (T+N) with respect to theta1
DTNQ_wrt_theta1 = BuildDTQ(fww_wrt_th1,alpha(ite).*gww_wrt_th1,k1,k2);

% Create The matrix (T+N) with respect to theta2
DTNQ_wrt_theta2 = BuildDTQ(fww_wrt_th2,alpha(ite).*gww_wrt_th2,k1,k2);

%%
% Create the matrix C for input into iteration

H_z     = DYG-DPG;

H_x     = DTNQ*M;

H_alpha  = DTNQ_wrt_alpha*M*x_ls - ...
    (ck_wrt_alpha + ht_wrt_alpha);

H_theta1 = DTNQ_wrt_theta1*M*x_ls - ...
    (ck_wrt_theta1 + ht_wrt_theta1);

H_theta2 = DTNQ_wrt_theta2*M*x_ls - ...
    (ck_wrt_theta2 + ht_wrt_theta2);

C       = [H_z H_x H_alpha H_theta1 H_theta2];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    zk;...
    x_ls;...
    alpha(ite);...
    th1(ite);...
    th2(ite)
    ];

yy              =   start_point;

% Set the termination criterion to a large value. It will be
% over written later.
condition(ite) = norm(res_vec)/norm(ck);

% Edit 17/11/2015
xk = x_ls;

while (condition(ite) >(SETTINGS.MAX_ERROR_SNTLN) &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN)
    
    % Use the QR decomposition to solve the LSE problem
    % min |y-p| subject to Cy=q
    
    y = LSE(E,p,C,res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % obtain the small changes
    
    % % Break down y into its sections
    
    % Get the entries corresponding to perturbations of f(x,y) and g(x,y)
    delta_zk        = y(1:nCoeff_fxy + nCoeff_gxy ,1);
    
    % Remove the entries in vector y corresponding to perturbations of f(x,y) and g(x,y).
    y(1:nCoeff_fxy + nCoeff_gxy) = [];
    
    % Remove the entries from the vector y
    delta_xk        = y(1:nCoeff_x,1);
    
    % Remove them from the list of coefficients
    y(1:nCoeff_x) = [];
    
    % Get the entry corresponding to alpha
    delta_alpha     = y(1:1);
    
    % Remove the entry from vector y
    y(1) = [];
    
    % Get the entry corresponding to theta1
    delta_th1    = y(1:1);
    
    % Remove the entry from vector y
    y(1) = [];
    
    % Get the coefficient corresponding to theta2
    delta_theta2    = y(1:1);

    % % Update the variables
    
    % Update variables z_{k}, where z_{k} are perturbations in the
    % coefficients of f and g.
    zk = zk + delta_zk;
    
    % Update x_{k}, where x_{k} is the solution vector, containing
    % coefficients u and v.
    xk = xk + delta_xk;
    
    % Update alpha
    alpha(ite) = alpha(ite-1) + delta_alpha;
    
    % Update theta_{1}
    th1(ite) = th1(ite-1) + delta_th1;
    
    % Update theta_{2}
    th2(ite) = th2(ite-1) + delta_theta2;
    
    % %
    % Obtain polynomials in modified bersntein basis a_{i}\theta^{i}
    
    % Obtain new f(w,w) with improved theta1, and theta2
    fww_matrix = GetWithThetas(fxy_matrix,th1(ite),th2(ite));
    
    % Obtain new g(w,w) with improved theta1 and theta2
    gww_matrix = GetWithThetas(gxy_matrix,th1(ite),th2(ite));
        
    % Construct the Sylvester subresultant matrix DTQ.
    DTQ_fg = BuildDTQ(fww_matrix, alpha(ite).*gww_matrix,k1,k2);
    
    % %
    % Get partial derivatives
    
    % Get the partial derivative of f(\omega_{1},\omega_{2}) with respect 
    % to \alpha
    fww_wrt_alpha            = zeros(m1+1,m2+1);
    
    % Get the partial derivative of g(\omega_{1},\omega_{2}) with respect
    % to \alpha
    alpha_gww_wrt_alpha      = gww_matrix;
    
    % Calculate the partial derivatives of fw and gw with respect to theta1
    % divide the rows by theta1 and multiply by the old power
    
    % Get the partial derivative of f(w1,w2) with respect to theta_{1}
    fww_wrt_th1 = Differentiate_wrt_theta1(fww_matrix,th1(ite));
    
    % Get the partial derivative of g(w1,w2) with respect to theta_{1}
    gww_wrt_th1 = Differentiate_wrt_theta1(gww_matrix,th1(ite));
        
    % Get the partial derivative of f with respect to theta_{2}
    fww_wrt_th2 = Differentiate_wrt_theta2(fww_matrix,th2(ite));
    
    % Get the partial derivative of g with respect to theta_{2}
    gww_wrt_th2 = Differentiate_wrt_theta2(gww_matrix,th2(ite));
    
    % Calculate the Partial derivative of T with respect to alpha.
    DTQ_wrt_alpha = BuildDTQ(fww_wrt_alpha, alpha_gww_wrt_alpha,k1,k2);
        
    % Calculate the partial derivative of DTQ with respect to theta_{1}
    DTQ_wrt_th1 = BuildDTQ(fww_wrt_th1,alpha(ite).*gww_wrt_th1,k1,k2);
    
    % Calculate the partial derivative of DTQ with respect to theta_{2}
    DTQ_wrt_th2 = BuildDTQ(fww_wrt_th2,alpha(ite).*gww_wrt_th2,k1,k2);
    
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ck = DTQ_fg*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha, \theta_{1}
    % and \theta_{2}
    ck_wrt_alpha    = DTQ_wrt_alpha*e;
    ck_wrt_theta1   = DTQ_wrt_th1*e;
    ck_wrt_theta2   = DTQ_wrt_th2*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fx      = zk(1:nCoeff_fxy);
    z_gx      = zk(nCoeff_fxy + 1 :end);
    
    % Get the vectors z_fx and z_gx as matrices, which match the shape of
    % f(x) and g(x).
    z_fx_mat = GetAsMatrix(z_fx,m1,m2);
    z_gx_mat = GetAsMatrix(z_gx,n1,n2);
    
    % Get matrices z_fw_mat and z_gw_mat, by multiplying rows by
    % theta_{1}^{i} and columns by theta_{2}^{j}
    z_fw_mat = GetWithThetas(z_fx_mat,th1(ite),th2(ite));
    z_gw_mat = GetWithThetas(z_gx_mat,th1(ite),th2(ite));
    
    % Calculate the derivatives of z_fw and z_gw with repect to \alpha.
    zfw_wrt_alpha    = zeros(m1+1,m2+1);
    alpha_zgw_wrt_alpha    = z_gw_mat;
    
    % Calculate the derivative of z_fw with respect to \theta_{1}.
    zfw_wrt_theta1 = Differentiate_wrt_theta1(z_fw_mat,th1(ite));
    
    % Calculate the derivative of z_fw with respect to \theta_{2}
    zfw_wrt_theta2 = Differentiate_wrt_theta2(z_fw_mat,th2(ite));
    
    % Calculate the derivative of z_gw with respect ot theta1
    zgw_wrt_theta1 = Differentiate_wrt_theta1(z_gw_mat,th1(ite));
    
    % Calculate the deriviate of z_gw with respect to theta2
    zgw_wrt_theta2 = Differentiate_wrt_theta2(z_gw_mat,th2(ite));
    
    % Build the coefficient Matrix N = [T(z1) T(z2)], of structured perturbations, with
    % same structure as DTQ.
    DNQ = BuildDTQ(z_fw_mat,alpha(ite).*z_gw_mat,k1,k2);
    
    % Build the coefficient matrix N with respect to alpha
    DNQ_wrt_alpha = BuildDTQ(zfw_wrt_alpha, alpha_zgw_wrt_alpha,k1,k2);
    
    % Calculate the derivatives of DNQ with respect to theta
    DNQ_wrt_theta1 = BuildDTQ(zfw_wrt_theta1, alpha(ite).*zgw_wrt_theta1,k1,k2);
    
    % Calculate the derivatives of DNQ with respect to theta
    DNQ_wrt_theta2 = BuildDTQ(zfw_wrt_theta2, alpha(ite).*zgw_wrt_theta2,k1,k2);
        
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{k} the column of S_{k} moved to the RHS
    ht = DNQ*e;
    
    % Calculate the derivative of h with respect to alpha
    h_alpha = DNQ_wrt_alpha*e;
    
    % Calculate the derivative of h with respect to theta1
    h_theta1 = DNQ_wrt_theta1*e;
    
    % Calculate the derivative of h with respect to theta2
    h_theta2 = DNQ_wrt_theta2*e;
    
    % Build the matrix (T+N)
    DTNQ = BuildDTQ(...
                    fww_matrix + z_fw_mat,...
                    alpha(ite).*(gww_matrix + z_gw_mat),...
                    k1,k2);
    
    % Calculate the paritial derivative of (T+N) with respect to
    % alpha
    DTNQ_alpha = BuildDTQ(...
                    fww_wrt_alpha + zfw_wrt_alpha,...
                    alpha_gww_wrt_alpha + alpha_zgw_wrt_alpha,...
                    k1,k2);
    
    
    % Calculate the paritial derivative of (T+N) with respect to theta1
    DTNQ_theta1 = BuildDTQ(...
                    fww_wrt_th1 + zfw_wrt_theta1,...
                    alpha(ite).*gww_wrt_th1 + zgw_wrt_theta1,...
                    k1,k2);
    
    % Calculate the paritial derivative of (T+N) with respect to theta2
    DTNQ_theta2 = BuildDTQ(...
                    fww_wrt_th2 + zfw_wrt_theta2,...
                    alpha(ite).*gww_wrt_th2 + zgw_wrt_theta2,...
                    k1,k2);
    
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    Y = BuildY(m1,m2,n1,n2,k1,k2,idx_col,xk,alpha(ite),th1(ite),th2(ite));
    % Build DYG
    DYG = D*Y*G;
       
    % Calculate the matrix DP where P is the matrix such that c = P[f;g]
    P = BuildP(m1,m2,n1,n2,th1(ite),th2(ite),alpha(ite),k1,k2,idx_col);
    % Build DPG
    DPG = D*P*G;
    
    % Get residual as a vector
    rk = (ck+ht) - DTNQ*M*xk ;
    
    % Create the matrix C. This is made up of five submatrices, HZ, Hx,
    % H_alpha and H_theta1 and H_theta2.
    
    Hz          = DYG - DPG;
    
    Hx          = DTNQ*M;
    
    H_alpha     = DTNQ_alpha*M*xk - (ck_wrt_alpha + h_alpha);
    
    H_theta1    = DTNQ_theta1*M*xk - (ck_wrt_theta1 + h_theta1);
    
    H_theta2    = DTNQ_theta2*M*xk - (ck_wrt_theta2 + h_theta2);
    
    C = [Hz,Hx,H_alpha,H_theta1, H_theta2];  % the matrix C
    
    % Calculate the new right hand vector
    ek = ck + ht;
    
    % Update gnew - used in lse problem
    res_vec = rk;
    
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(rk) / norm(ek);
    
    % Update fnew - used in LSE Problem.
    p = -(yy-start_point);
    
    
end

% Plot Graphs

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure('name','SNTLN - Residuals in SNTLN')
        hold on
        title('Residuals in SNTLN')
        xlabel('Iterations')
        ylabel('log_{10} Residuals')
        plot((1:1:ite),log10(condition),'-s')
        hold off

    case 'n'
        % Do Nothing
    otherwise
        error('err')
end


%

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.

% get the vector zk

zPert_f_vec = zk(1:nCoeff_fxy);
zPert_f_mat = GetAsMatrix(zPert_f_vec,m1,m2);

zPert_g_vec = zk(nCoeff_fxy+1:end);
zPert_g_mat = GetAsMatrix(zPert_g_vec,n1,n2);

fxy_output = fxy_matrix_n + zPert_f_mat;

gxy_output = gxy_matrix_n + zPert_g_mat;

X_output  = xk;

alpha_output = alpha(ite);

theta1_output = th1(ite);

theta2_output = th2(ite);


% Print the number of iterations
fprintf('--------------------------------------------------------------------------- \n')
fprintf('Iterations over SNTLN fucnction for low rank approximation : %i \n', ite);
fprintf('--------------------------------------------------------------------------- \n')
end

















