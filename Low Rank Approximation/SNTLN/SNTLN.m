function [ fxy_lr,gxy_lr,uxy_lr,vxy_lr,alpha_lr,th1_lr,th2_lr] = ...
    SNTLN( fxy,gxy, i_alpha, i_th1, i_th2,k1,k2, idx_col)
% Obtain the low rank approximation by method of SNTLN (Structured
% Non-Linear Total Least Norm) of the Sylvester matrix S(f,g).
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) in the Bernstein basis.
%
% gxy : (Matrix) Coefficients of polynomial g(x,y) in the Bernstein basis.
%
% i_alpha : (Float) Initial value of alpha
%
% i_theta : (Float) Initial value of theta
%
% k1 : (Int) Degree of polynomial d(x,y) with respect to x
%
% k2 : (Int) Degree of polynomial d(x,y) with respect to y
%
% idx_col : (Int) Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.
%
% % Outputs.
%
% fxy_lr : (Matrix) Coefficients of f(x,y) in the Bernstein basis,
% perturbations.
%
% gxy_lr : (Matrix) Coefficients of g(x,y) in the Bernstein basis with added 
% perturbations.
%
% uxy_lr : (Matrix) 
%
% vxy_lr : (Matrix) 
%
% alpha_lr : (Float) 
%
% th1_lr : (Float) 
%
% th2_lr : (Float)



% Global Inputs
global SETTINGS

if isempty(SETTINGS.MAX_ERROR_SNTLN) || isempty(SETTINGS.MAX_ITERATIONS_SNTLN)
    error('err')
end


% Set the initial iterations number
ite = 1;

% Set initial values of alpha and theta
th1(ite) = i_th1;
th2(ite) = i_th2;
alpha(ite) = i_alpha;

% Get degree of polynomials f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of polynomial g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Get the number of coefficients in the polynomial f(x,y)
nCoefficients_fxy = (m1 + 1) * (m2 + 1);

% Get the number of coefficients in the polynomial g(x,y)
nCoefficients_gxy = (n1 + 1) * (n2 + 1);

% Get the number of coefficients in both f(x,y) and g(x,y)
nCoefficients_fg = nCoefficients_fxy + nCoefficients_gxy;

% Get the number of coefficients in v(x,y)
nCoefficients_vxy = (n1 - k1 + 1) * (n2 - k2 + 1);

% Get the number of coefficients in u(x,y)
nCoefficients_uxy = (m1 - k1 + 1) * (m2 - k2 + 1);

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
nCoefficients_x = nCoefficients_uxy + nCoefficients_vxy - 1;

% Get the number of columns in C_{n-k}(f), the first partition of the Sylvester
% Matrix S_{t}(f,g).
nColumns_Tf = nCoefficients_vxy;

% Get the number of columns in C_{t}(g), the second partition of the
% Sylvester matrix S_{t}(f,g)
nColumns_Tg = nCoefficients_uxy;

% Get the total number of columns in the Sylvester matrix S_{t}(f,g)
nColumns_Sk = nColumns_Tf + nColumns_Tg;

% Create the identity matrix
I = eye(nColumns_Sk, nColumns_Sk);

% Create the matrix M, such that S(f,g)*M gives A_{t}, the Sylvester Matrix
% with the optimal column removed.
M = I;
M(:, idx_col) = [];

% Let e be the column removed from the identity matrix, such that
% S_{t}(f,g) * e gives the column c_{t}, where c_{t} is the optimal column
% removed from the Sylvester subresultant.
e = I(:, idx_col);

% % Preprocessing

% Get the polynomial f(\omega_{1},\omega_{2})
fww = GetWithThetas(fxy, th1, th2);

% Get the polynomial g(\omega_{1},\omega_{2})
gww = GetWithThetas(gxy, th1, th2);

% % Form the Coefficient Matrix T = [C(f)|C(g)] such that DTQ * x = [col]
DTQ_fg = BuildDTQ_Bivariate_2Polys(fww, alpha.*gww, k1, k2);

% Get partial derivative of f(\omega_{1},\omega_{2}) with respect to \alpha
fww_wrt_alpha            = zeros(m1+1, m2+1);

% Get partial derivative of g(\omega_{1},\omega_{2}) with respect to \alpha
alpha_gww_wrt_alpha      = gxy;

% Calculate the partial derivatives of f(w,w) with respect to theta_1
fww_wrt_th1 = Differentiate_wrt_th1(fww, th1(ite));

% Calculate the partial derivative of g(w,w) with respect to theta_1
gww_wrt_th1 = Differentiate_wrt_th1(gww, th1(ite));

% Calculate the partial derivative of f(w,w) with respect to theta_2
fww_wrt_th2 = Differentiate_wrt_th2(fww, th2(ite));

% Calculate the partial deriviates of g(w,w) with respect to theta_2
gww_wrt_th2 = Differentiate_wrt_th2(gww, th2(ite));

% Build the derivative of T(f,g) with respect to alpha
DTQ_alpha = BuildDTQ_Bivariate_2Polys(fww_wrt_alpha, alpha_gww_wrt_alpha, k1, k2);

% Calculate the derivative of T(f,g) with respect to theta_{1}
DTQ_wrt_th1 = BuildDTQ_Bivariate_2Polys(fww_wrt_th1, alpha(ite).*gww_wrt_th1, k1, k2);

% Calcualte the derivative of T(f,g) with respect to theta_2
DTQ_wrt_th2 = BuildDTQ_Bivariate_2Polys(fww_wrt_th2, alpha(ite).*gww_wrt_th2, k1, k2);


%
% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.

nRows_Sylvester_matrix = (m1 + n1 - k1 + 1) * (m2 + n2 - k2 + 1);

zk = zeros(nCoefficients_fg , 1);

%
% Initilaise the derivative of N wrt alpha.
DNQ_wrt_alpha = zeros(nRows_Sylvester_matrix, nColumns_Tf + nColumns_Tg);

% Initilaise the derivative of N wrt theta_1.
DNQ_wrt_th1 = zeros(nRows_Sylvester_matrix, nColumns_Tf + nColumns_Tg);

% Initialise the derivative of N wrt theta 2
DNQ_wrt_th2 = zeros(nRows_Sylvester_matrix, nColumns_Tf + nColumns_Tg);

%
% Initialise the derivative of h
% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
hk_wrt_alpha = DNQ_wrt_alpha*e;
hk_wrt_th1 = DNQ_wrt_th1*e;
hk_wrt_th2 = DNQ_wrt_th2*e;

% Get the matrix A_{k}(f,g), which is the subresultant matrix S(f,g) with
% an opitmal column removed
Ak_fg = DTQ_fg;
ck = DTQ_fg(:, idx_col);
Ak_fg(:, idx_col) = [];

% % Build the matrix D*P*G

DPG = BuildDPG_SNTLN(m1, m2, n1, n2, th1(ite), th2(ite), alpha(ite), k1, k2, idx_col);

% Test
f = GetAsVector(fxy);
g = GetAsVector(gxy);
test1a = DPG * [f;g];
test1b = ck;
display(norm(test1a - test1b));

%
% Calculate the derivatives wrt alpha and theta of the removed column.
ck_wrt_alpha = DTQ_alpha * e;
ck_wrt_th1 = DTQ_wrt_th1 * e;
ck_wrt_th2 = DTQ_wrt_th2 * e;


%
% Perform QR decomposition of Ak to obtain the solution x
xk = SolveAx_b(Ak_fg, ck);

% Insert a zero into the position of the optimal_column.
% Partition x_ls into the two parts for coefficients of 

first_part = xk( 1 : (idx_col - 1) );
second_part = xk( idx_col : end );
x = [first_part ; 0 ; second_part];

% Build Matrix DYG, where D*Y(x1,x2)*G*[f;g] = D*T(f,g)*Q*[u;v]
DYG = BuildDYG_SNTLN(m1, m2, n1, n2, k1, k2, x, alpha(ite), th1(ite), th2(ite));


% Test
test2a = DTQ_fg * x;
test2b = DYG * [f;g];
test2 = test2a - test2b;
display(norm(test2))

% Calculate the initial residual r = ck - (Ak*x)
vResidual = ck - (DTQ_fg * M * xk);

% % Get the matrix p, which will store all the perturbations returned 
% from LSE file
nEntries = nCoefficients_fxy...
    + nCoefficients_gxy ...
    + nCoefficients_x ...
    + 3;


% Set the intial value of E to the identity matrix
E = eye(nEntries);

% Create the matrix D(T+N)Q, initially N is empty so this is the same as T.
DTNQ = BuildDTQ_Bivariate_2Polys(fww, alpha(ite).*gww, k1, k2);

% Create The matrix (T+N) with respect to alpha
DTNQ_wrt_alpha = BuildDTQ_Bivariate_2Polys(fww_wrt_alpha, alpha_gww_wrt_alpha, k1, k2);

% Create The matrix (T+N) with respect to theta1
DTNQ_wrt_theta1 = BuildDTQ_Bivariate_2Polys(fww_wrt_th1, alpha(ite).*gww_wrt_th1, k1, k2);

% Create The matrix (T+N) with respect to theta2
DTNQ_wrt_theta2 = BuildDTQ_Bivariate_2Polys(fww_wrt_th2, alpha(ite).*gww_wrt_th2, k1, k2);

% %
% Create the matrix C for input into iteration

H_z     = DYG - DPG;

H_x     = DTNQ*M;

H_alpha  = DTNQ_wrt_alpha*M*xk - ...
    (ck_wrt_alpha + hk_wrt_alpha);

H_th1 = DTNQ_wrt_theta1*M*xk - ...
    (ck_wrt_th1 + hk_wrt_th1);

H_th2 = DTNQ_wrt_theta2*M*xk - ...
    (ck_wrt_th2 + hk_wrt_th2);

C = [H_z H_x H_alpha H_th1 H_th2];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
        zk;...
        xk;...
        alpha(ite);...
        th1(ite);...
        th2(ite)
    ];

yy  = start_point;

f = -(yy - start_point);

% Set the termination criterion to a large value. It will be
% over written later.
condition(ite) = norm(vResidual) / norm(ck);

while (condition(ite) >(SETTINGS.MAX_ERROR_SNTLN) &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN)
    
    % Use the QR decomposition to solve the LSE problem
    % min |y-p| subject to Cy=q
    
    y = LSE(E, f, C, vResidual);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % Get the entries corresponding to perturbations of f(x,y) and g(x,y)
    delta_zk = y(1 : nCoefficients_fxy + nCoefficients_gxy ,1);
    
    % Remove the entries in vector y corresponding to perturbations of f(x,y) and g(x,y).
    y(1 : nCoefficients_fxy + nCoefficients_gxy) = [];
    
    % Get the coefficients of xk
    delta_xk        = y(1 : nCoefficients_x, 1);
    
    % Remove them from the list of coefficients
    y(1 : nCoefficients_x) = [];
    
    % Get the entry corresponding to alpha
    delta_alpha     = y(1 : 1);
    
    % Remove the entry from vector y
    y(1) = [];
    
    % Get the entry corresponding to theta1
    delta_th1    = y(1:1);
    
    % Remove the entry from vector y
    y(1) = [];
    
    % Get the coefficient corresponding to theta2
    delta_th2    = y(1:1);

    % % Update the variables
    
    % Update variables z_{k}, where z_{k} are perturbations in the
    % coefficients of f and g.
    zk = zk + delta_zk;
    
    % Update x_{k}, where x_{k} is the solution vector, containing
    % coefficients u and v.
    xk = xk + delta_xk;
    
    % Update alpha
    alpha(ite) = alpha(ite - 1) + delta_alpha;
    
    % Update theta_{1}
    th1(ite) = th1(ite - 1) + delta_th1;
    
    % Update theta_{2}
    th2(ite) = th2(ite - 1) + delta_th2;
    
    % %
    % Obtain polynomials in modified bersntein basis a_{i}\theta^{i}
    
    % Obtain new f(\omega_{1},\omega_{2}) and g(\omega_{1},\omega_{2} with 
    % improved \theta_{1}, and \theta_{2}
    fww = GetWithThetas(fxy, th1(ite), th2(ite));
    gww = GetWithThetas(gxy, th1(ite), th2(ite));
        
    % Construct the kth Sylvester subresultant matrix DTQ.
    DTQ_fg = BuildDTQ_Bivariate_2Polys(fww, alpha(ite).*gww, k1, k2);
    
    % %
    % Get partial derivatives
    
    % Get the partial derivative of f(\omega_{1},\omega_{2}) with respect 
    % to \alpha
    fww_wrt_alpha            = zeros(m1 + 1, m2 + 1);
    
    % Get the partial derivative of g(\omega_{1},\omega_{2}) with respect
    % to \alpha
    alpha_gww_wrt_alpha      = gww;
    
    % Calculate the partial derivatives of fw and gw with respect to theta1
    % divide the rows by theta1 and multiply by the old power
    
    % Get the partial derivative of f and g with respect to theta_{1}
    fww_wrt_th1 = Differentiate_wrt_th1(fww, th1(ite));
    gww_wrt_th1 = Differentiate_wrt_th1(gww, th1(ite));
        
    % Get the partial derivative of f and g with respect to theta_{2}
    fww_wrt_th2 = Differentiate_wrt_th2(fww, th2(ite));
    gww_wrt_th2 = Differentiate_wrt_th2(gww, th2(ite));
    
    % Calculate the Partial derivative of T with respect to alpha.
    DTQ_wrt_alpha = BuildDTQ_Bivariate_2Polys(fww_wrt_alpha, alpha_gww_wrt_alpha, k1, k2);
        
    % Calculate the partial derivative of DTQ with respect to theta_{1}
    DTQ_wrt_th1 = BuildDTQ_Bivariate_2Polys(fww_wrt_th1, alpha(ite).*gww_wrt_th1, k1, k2);
    
    % Calculate the partial derivative of DTQ with respect to theta_{2}
    DTQ_wrt_th2 = BuildDTQ_Bivariate_2Polys(fww_wrt_th2, alpha(ite).*gww_wrt_th2, k1, k2);
    
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ck = DTQ_fg*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha, \theta_{1}
    % and \theta_{2}
    ck_wrt_alpha = DTQ_wrt_alpha*e;
    ck_wrt_th1   = DTQ_wrt_th1*e;
    ck_wrt_th2   = DTQ_wrt_th2*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    vec_z_fx      = zk(1 : nCoefficients_fxy);
    vec_z_gx      = zk(nCoefficients_fxy + 1 : end);
    
    % Get the vectors z_fx and z_gx as matrices, which match the shape of
    % f(x) and g(x).
    z_fxy = GetAsMatrix(vec_z_fx, m1, m2);
    z_gxy = GetAsMatrix(vec_z_gx, n1, n2);
    
    % Get matrices z_fw_mat and z_gw_mat, by multiplying rows by
    % theta_{1}^{i} and columns by theta_{2}^{j}
    z_fww = GetWithThetas(z_fxy, th1(ite), th2(ite));
    z_gww = GetWithThetas(z_gxy, th1(ite), th2(ite));
    
    % Calculate the derivatives of z_fw and z_gw with repect to \alpha.
    zfw_wrt_alpha           = zeros(m1+1,m2+1);
    alpha_zgw_wrt_alpha     = z_gww;
    
    % Calculate the derivative of z_fw and z_gw with respect to \theta_{1}.
    zfw_wrt_th1 = Differentiate_wrt_th1(z_fww, th1(ite));
    zgw_wrt_th1 = Differentiate_wrt_th1(z_gww, th1(ite));
    
    % Calculate the derivative of z_fw and z_gw with respect ot \theta_{1}
    zgw_wrt_th2 = Differentiate_wrt_th2(z_gww, th2(ite));
    zfw_wrt_th2 = Differentiate_wrt_th2(z_fww, th2(ite));
    
    % Build the coefficient Matrix N = [T(z1) T(z2)], of structured perturbations, with
    % same structure as DTQ.
    DNQ = BuildDTQ_Bivariate_2Polys(z_fww, alpha(ite).*z_gww, k1, k2);
    
    % Build the coefficient matrix N with respect to alpha
    DNQ_wrt_alpha = BuildDTQ_Bivariate_2Polys(zfw_wrt_alpha, alpha_zgw_wrt_alpha, k1, k2);
    
    % Calculate the derivatives of DNQ with respect to theta
    DNQ_wrt_th1 = BuildDTQ_Bivariate_2Polys(zfw_wrt_th1, alpha(ite).*zgw_wrt_th1, k1, k2);
    
    % Calculate the derivatives of DNQ with respect to theta
    DNQ_wrt_th2 = BuildDTQ_Bivariate_2Polys(zfw_wrt_th2, alpha(ite).*zgw_wrt_th2, k1, k2);
        
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{k} the column of S_{k} moved to the RHS
    hk = DNQ * e;
     
    % Calculate the derivative of h with respect to alpha
    h_alpha = DNQ_wrt_alpha * e;
    
    % Calculate the derivative of h with respect to theta1
    h_th1 = DNQ_wrt_th1 * e;
    
    % Calculate the derivative of h with respect to theta2
    h_th2 = DNQ_wrt_th2 * e;
    
    % Build the matrix (T+N)
    DTNQ = BuildDTQ_Bivariate_2Polys(...
                    fww + z_fww,...
                    alpha(ite).*(gww + z_gww),...
                    k1, k2);
    
    % Calculate the paritial derivative of (T+N) with respect to
    % alpha
    DTNQ_alpha = BuildDTQ_Bivariate_2Polys(...
                    fww_wrt_alpha + zfw_wrt_alpha,...
                    alpha_gww_wrt_alpha + alpha_zgw_wrt_alpha,...
                    k1, k2);
    
    
    % Calculate the paritial derivative of (T+N) with respect to theta1
    DTNQ_th1 = BuildDTQ_Bivariate_2Polys(...
                    fww_wrt_th1 + zfw_wrt_th1,...
                    alpha(ite).*gww_wrt_th1 + zgw_wrt_th1,...
                    k1, k2);
    
    % Calculate the paritial derivative of (T+N) with respect to theta2
    DTNQ_th2 = BuildDTQ_Bivariate_2Polys(...
                    fww_wrt_th2 + zfw_wrt_th2,...
                    alpha(ite).*gww_wrt_th2 + zgw_wrt_th2,...
                    k1, k2);
    
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    % Insert a zero into the position of the optimal_column.
    % Partition x_ls into the two parts for coefficients of 
    first_part = xk(1 : (idx_col-1));
    second_part = xk(idx_col : end);
    x = [first_part ; 0 ; second_part];
    
    DYG = BuildDYG_SNTLN(m1, m2, n1, n2, k1, k2, x, alpha(ite), th1(ite), th2(ite));
    
    
    % Calculate the matrix DP where P is the matrix such that c = P[f;g]
    DPG = BuildDPG_SNTLN(m1, m2, n1, n2, th1(ite), th2(ite), alpha(ite), k1, k2, idx_col);
       
    
    % Get residual as a vector
    vResidual = (ck + hk) - DTNQ*M*xk ;
    
    % Create the matrix C. This is made up of five submatrices, HZ, Hx,
    % H_alpha and H_theta1 and H_theta2.
    
    Hz = DYG - DPG;
    
    Hx = DTNQ*M;
    
    H_alpha = DTNQ_alpha*M*xk - (ck_wrt_alpha + h_alpha);
    
    H_th1 = DTNQ_th1*M*xk - (ck_wrt_th1 + h_th1);
    
    H_th2 = DTNQ_th2*M*xk - (ck_wrt_th2 + h_th2);
    
    C = [Hz, Hx, H_alpha, H_th1, H_th2];  % the matrix C
    
        
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(vResidual) / norm(ck + hk);
    
    % Update fnew - used in LSE Problem.
    f = -(yy - start_point);
    
    
end

% Print the number of iterations
LineBreakLarge()
fprintf([mfilename ' : ' sprintf('Iterations over SNTLN fucnction for low rank approximation : %i \n', ite)]);
LineBreakLarge()
SETTINGS.LOW_RANK_APPROX_REQ_ITE = ite;


% Plot Graphs
PlotSNTLN();


%

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.

% get the vector zk
zPert_f_vec = zk(1 : nCoefficients_fxy);
zPert_f_mat = GetAsMatrix(zPert_f_vec,m1,m2);

zPert_g_vec = zk(nCoefficients_fxy+1 : end);
zPert_g_mat = GetAsMatrix(zPert_g_vec, n1, n2);

% Get f(x,y) low rank
fxy_lr = fxy + zPert_f_mat;

% Get g(x,y) low rank
gxy_lr = gxy + zPert_g_mat;

% %
% Get u(x,y) and v(x,y)
first_part = xk(1:(idx_col-1));
second_part = xk(idx_col:end);
x = [first_part ; -1 ; second_part];

% Get vectors of u(\omega_{1},\omega_{2})
v_vww = x(1:nCoefficients_vxy);
v_uww = -1.*x(nCoefficients_vxy+1:end);

% Get v(\omega_{1},\omega_{2}) and u(\omega_{1},\omega_{2})
vww = GetAsMatrix(v_vww, n1-k1, n2-k2);
uww = GetAsMatrix(v_uww, m1-k1, m2-k2);

% Get v(x,y) and u(x,y)
vxy_lr = GetWithoutThetas(vww, th1(ite), th2(ite));
uxy_lr = GetWithoutThetas(uww, th1(ite), th2(ite));

% %
% Get \alpha, \theta_{1}
alpha_lr = alpha(ite);

th1_lr = th1(ite);

th2_lr = th2(ite);


end

















