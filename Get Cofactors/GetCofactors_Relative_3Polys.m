function [uxy, vxy, wxy] = GetCofactors_Relative_3Polys(fxy, gxy, hxy, t1, t2)
% Given two input polynomials and the degree of the GCD return the quotient
% polynomials u(x,y) and v(x,y)
%
% % Inputs
%
% [fxy, gxy, wxy] : Coefficients of polynomial f(x,y), g(x,y) and h(x,y)
%
% t1 : Degree of d(x,y) with respect to x
%
% t2 : Degree of d(x,y) with respect to y
%
% % Outputs
%
% [uxy, vxy, wxy] : Coefficients of polynomial u(x,y), v(x,y) and w(x,y)


% Initialise Global Variables
global SETTINGS

% % 
%

% Get the degrees of polynomial f(x,y)
[m1, m2] = GetDegree(fxy);

% Get the degrees of polynomial g(x,y)
[n1, n2] = GetDegree(gxy);

% Get the degrees of polynomial h(x,y)
[o1, o2] = GetDegree(hxy);

% Build the (t1,t2)-th subresultant
Sk1k2 = BuildDTQ_3Polys(fxy, gxy, hxy, t1, t2);


%% Find Optimal column for removal from St
% given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized
opt_col_index = GetOptimalColumn(Sk1k2);

% Get the matrix A_{t_{1},t_{2}} 
At = Sk1k2;
At(:,opt_col_index) = [];

% Get the vector c_{t_{1},t_{2}} removed from S_{t_{1},t_{2}}
ct = Sk1k2(:,opt_col_index);


%% Get the coefficients for u(x,y) and v(x,y)
x_ls = SolveAx_b(At,ct);


% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col_index)-1);
    -1;
    x_ls(opt_col_index:end);
    ];

nCoeffs_vxy = (n1-t1+1) * (n2-t2+1);
nCoeffs_uxy = (m1-t1+1) * (m2-t2+1);
nCoeffs_wxy = (o1-t1+1) * (o2-t2+1);


% Get coefficients of u(x,y) as a vector 
v_vxy = vecx(1:nCoeffs_vxy);

v_wxy = vecx(nCoeffs_vxy + 1 : nCoeffs_vxy + nCoeffs_wxy];

% Get coefficients of v(w,w) as a vector
v_uxy = -1 .* vecx(nCoeffs_vxy + nCoeffs_wxy + 1 : end);


% Get u(x,y) as a matrix of coefficients
uxy = GetAsMatrix(v_uxy, m1-t1, m2-t2);

% Get v(x,y) as a matrix of coefficients
vxy = GetAsMatrix(v_vxy, n1-t1, n2-t2);

% Get w(x,y) as a matrix of coefficients
wxy = GetAsMatrix(v_wxy, o1-t1, o2-t2);

switch SETTINGS.SYLVESTER_BUILD_METHOD
    case 'T'
        
        vxy = GetWithoutBinomails(vxy);
        uxy = GetWithoutBinomials(uxy);
        wxy = GetWithoutBinomials(wxy);
        
    case 'DT'
        
        vxy = GetWithoutBinomails(vxy);
        uxy = GetWithoutBinomials(uxy);
        wxy = GetWithoutBinomials(wxy);
        
    case 'DTQ'
        
    case 'TQ'
        
        

    otherwise
        error('err')
end

end