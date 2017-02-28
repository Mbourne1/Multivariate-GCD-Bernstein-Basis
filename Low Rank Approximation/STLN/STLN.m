function [ fxy_lr,gxy_lr,uxy_lr,vxy_lr] = STLN(fxy, gxy, k1, k2, idx_col)
% Obtain the low rank approximation by method of SNTLN (Structured
% Non-Linear Total Least Norm) of the Sylvester matrix S(f,g).
%
% % Inputs.
%
% fxy : Coefficients of polynomial f(x,y), in the Bernstein basis.
%
% gxy : Coefficients of polynomial g(x,y), in the Bernstein basis.
%
% k1 : Degree of AGCD with respect to x.
%
% k2 : Degree of AGCD with respect to y.
%
% idx_col : Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.
%
% % Outputs.
%
% fxy_lr : Coefficients of f(x,y) in the Bernstein basis, including added 
% perturbations.
%
% gxy_lr : Coefficients of g(x,y) in the Bernstein basis, including added 
% perturbations.
%
% uxy_lr : Coefficients of  
%
% vxy_lr :


% Global Variables
global SETTINGS
if isempty(SETTINGS.MAX_ERROR_SNTLN) || isempty(SETTINGS.MAX_ITERATIONS_SNTLN)
    error('err')
end

% Set the initial iterations number
ite = 1;

% Get degree of polynomial f(x,y).
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of polynomial g(x,y).
[n1, n2] = GetDegree_Bivariate(gxy);

% Get the number of coefficients in the polynomial f(x,y)
nCoeffs_fxy = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
nCoeffs_gxy = (n1+1) * (n2+1);

% Get the number of coefficients in both f(x,y) and g(x,y)
nCoeffs_fg = nCoeffs_fxy + nCoeffs_gxy;

% Get the number of coefficients in v(x,y)
nCoeffs_vxy = (n1-k1+1) * (n2-k2+1);

% Get the number of coefficients in u(x,y)
nCoeffs_uxy = (m1-k1+1) * (m2-k2+1);

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
nCoeffs_xls = nCoeffs_uxy + nCoeffs_vxy - 1;

% Create the identity matrix I, the matrix M formed from I by removing the
% column equivalent to the optimal column for removal from the Sylvester
% subresultant matrix, so that S_{t}(f,g)*M = A_{t}, where A_{t} is the
% Sylvester subresultant matrix with the column removed.

% Get the number of columns in C_{t}(f), the first partition of the Sylvester
% Matrix S_{t}(f,g)
nCols_Cf = nCoeffs_vxy;

% Get the number of columns in C_{t}(g), the second partition of the
% Sylvester matrix S_{t}(f,g)
nCols_Cg = nCoeffs_uxy;

% Get the total number of columns in the Sylvester matrix S_{t}(f,g)
nCols_St1t2 = nCols_Cf + nCols_Cg;

% Create the identity matrix
I = eye(nCols_St1t2, nCols_St1t2);

% Create the matrix M, such that S(f,g)*M gives A_{t}, the Sylvester Matrix
% with the optimal column removed.
M = I;
M(:,idx_col) = [];

% Let e be the column removed from the identity matrix, such that
% S_{t}(f,g) * e gives the column c_{t}, where c_{t} is the optimal column
% removed from the Sylvester subresultant.
e = I(:,idx_col);

% Form the Coefficient Matrix T = [C(f)|C(g)] such that DTQ * x = [col]

% % Build the Sylvester matrix D*T(f,g)*Q
DTQ = BuildDTQ_Bivariate_2Polys(fxy, gxy, k1, k2);

% Get the matrix A_{t}(f,g), which is the subresultant matrix S(f,g) with
% an opitmal column removed
Ak = DTQ;
ck = DTQ(:,idx_col);
Ak(:,idx_col) = [];

% % 
% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.
zk = zeros(nCoeffs_fg , 1);

% %
% Build the matrix of binomials corresponding to polynomial f(x,y)
DPG = BuildDPG_SNTLN(m1, m2, n1, n2, 1, 1, 1, k1, k2, idx_col);

% Test
f = GetAsVector(fxy);
g = GetAsVector(gxy);
test1a = DPG * [f;g];
test1b = ck;
display(norm(test1a - test1b));

%
% Perform QR decomposition of Ak to obtain the solution x
xk = SolveAx_b(Ak,ck);

% % Build Matrix DYG, where D*Y(x1,x2)*G*[f;g] = D*T(f,g)*Q*[u;v]

% Insert a zero into the position of the optimal_column.
% Partition x_ls into the two parts for coefficients of 
first_part = xk(1:(idx_col-1));
second_part = xk(idx_col:end);
x = [first_part ; 0 ; second_part];
DYG = BuildDYG_SNTLN(m1,m2,n1,n2,k1,k2,x,1,1,1);

% Test
test2a = DTQ*x;
test2b = DYG*[f;g];
test2 = test2a - test2b;
display(norm(test2))

% %
% Calculate the initial residual r = ck - (Ak*x)
res_vec = ck - (DTQ*M*xk);

% Get the matrix p, which will store all the perturbations returned from LSE file
nEntries = nCoeffs_fxy...
    + nCoeffs_gxy ...
    + nCoeffs_xls ...
    + 3;


% Set the intial value of E to the identity matrix
E = eye(nEntries);

% %
% Create the matrix D(T+N)Q, initially N is empty so this is the same as T.
DTNQ = BuildDTQ_Bivariate_2Polys(fxy, gxy, k1, k2);

% %
% Create the matrix C for input into iteration

H_z     = DYG - DPG;

H_x     = DTNQ * M;

C       = [H_z H_x];

% Define the starting vector for the iterations for the LSE problem.
start_point = ...
    [...
        zk;...
        xk;...
    ];

yy = start_point;

f = -(yy-start_point);

% Set the termination criterion to a large value. It will be
% over written later.
condition(ite) = norm(res_vec)/norm(ck);

while condition(ite) >(SETTINGS.MAX_ERROR_SNTLN) &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
    % Use the QR decomposition to solve the LSE problem
    % min |y-p| subject to Cy=q
    
    y = LSE(E,f,C,res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
        
    % Get the coefficients of z_{f}(x,y) and z_{g}(x,y)
    delta_zk        = y(1:nCoeffs_fxy + nCoeffs_gxy ,1);
    
    % Update variables z_{k}, where z_{k} are perturbations in the
    % coefficients of f and g.
    zk = zk + delta_zk;
    
    % Get the vector of structured perturbations z_{f}(x,y) 
    vec_z_fxy      = zk(1:nCoeffs_fxy);
    
    % Get the vector of coefficients of z_{g}(x,y)
    vec_z_gxy      = zk(nCoeffs_fxy + 1 :end);
    
    % Get the vectors z_fx and z_gx as matrices, which match the shape of
    % f(x) and g(x).
    z_fxy = GetAsMatrix(vec_z_fxy,m1,m2);
    z_gxy = GetAsMatrix(vec_z_gxy,n1,n2);
       
    % Get the coefficients corresponding to x
    delta_xk        = y(nCoeffs_fg+1:end);
    
    % Update x_{k}, where x_{k}
    xk = xk + delta_xk;
    
    % Construct the Sylvester subresultant matrix DTQ.
    DTQ = BuildDTQ_Bivariate_2Polys(fxy, gxy, k1, k2);
       
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ck = DTQ*e;
    
    % Build the coefficient Matrix N = [T(z1) T(z2)], of structured perturbations, with
    % same structure as DTQ.
    DNQ = BuildDTQ_Bivariate_2Polys(z_fxy, z_gxy, k1, k2);

    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{k} the column of S_{k} moved to the RHS
    hk = DNQ*e;
       
    % Build the matrix (T+N)
    DTNQ = BuildDTQ_Bivariate_2Polys((fxy + z_fxy), (gxy + z_gxy), k1, k2);

    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    % Insert a zero into the position of the optimal_column.
    
    % Partition x_ls into the two parts for coefficients of 
    first_part = xk(1:(idx_col-1));
    second_part = xk(idx_col:end);
    x = [first_part ; 0 ; second_part];
   
    % Update the matrix Y and DYG
    DYG = BuildDYG_SNTLN(m1,m2,n1,n2,k1,k2,x,1,1,1);
  
    % Get residual as a vector
    res_vec = (ck + hk) - DTNQ*M*xk ;
    
    % Create the matrix C. This is made up of five submatrices, HZ, Hx,
    % H_alpha and H_theta1 and H_theta2.
    
    Hz          = DYG - DPG;
    
    Hx          = DTNQ*M;
    
    % Update the matrix C for LSE Problem
    C = [Hz Hx];  
    
    % Update the RHS vector 
    ek = ck + hk;
    
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(res_vec) / norm(ek);
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    
end


% Print the number of iterations
LineBreakLarge()
fprintf([mfilename ' : ' sprintf('Iterations over STLN fucnction for low rank approximation : %i \n', ite)]);
LineBreakLarge()
SETTINGS.LOW_RANK_APPROX_REQ_ITE = ite;

% Plot Graphs
PlotSNTLN()


% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.

% get the vector zk

zPert_f_vec = zk(1:nCoeffs_fxy);
zPert_f_mat = GetAsMatrix(zPert_f_vec,m1,m2);

zPert_g_vec = zk(nCoeffs_fxy+1:end);
zPert_g_mat = GetAsMatrix(zPert_g_vec,n1,n2);

% Get f(x,y) with added perturbations
fxy_lr = fxy + zPert_f_mat;

% Get g(x,y) with added perturbations
gxy_lr = gxy + zPert_g_mat;

% Get u(x,y) and v(x,y)
xa = xk(1:idx_col-1);
xb = xk(idx_col:end);
x = [xa ; -1 ; xb];

vec_vxy = x(1:nCoeffs_vxy);
vec_uxy = -1.*x(nCoeffs_vxy+1:end);

vxy_lr = GetAsMatrix(vec_vxy,n1-k1,n2-k2);
uxy_lr = GetAsMatrix(vec_uxy,m1-k1,m2-k2);


end

















