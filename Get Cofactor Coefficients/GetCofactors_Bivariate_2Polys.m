function [uxy, vxy] = GetCofactors_Bivariate_2Polys(fxy, gxy, t1, t2)
% Given two input polynomials and the degree of the GCD return the quotient
% polynomials u(x,y) and v(x,y)
%
% % Inputs
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% t1 : Degree of d(x,y) with respect to x
%
% t2 : Degree of d(x,y) with respect to y
%
% % Outputs
%
% uxy : Coefficients of polynomial u(x,y)
%
% vxy : Coefficients of polynomial v(x,y)

% Initialise Global Variables
global SETTINGS

% % 
%

% Get the degrees of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get the degrees of polynomial g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Build the (t1,t2)-th subresultant
Sk1k2 = BuildSubresultant_Bivariate_2Polys(fxy, gxy, t1, t2);

% Find Optimal column for removal from St
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

% Get number of coefficients in v(x,y) and u(x,y)
nCoeffs_vxy = (n1-t1+1) * (n2-t2+1);
nCoeffs_uxy = (m1-t1+1) * (m2-t2+1);


% Get coefficients of u(x,y) as a vector 
v_vxy = vecx(1:nCoeffs_vxy);

% Get coefficients of v(w,w) as a vector
v_uxy = -1 .* vecx(nCoeffs_vxy+1:end);

% % Get the value of fv-gu
norm(Sk1k2 * [v_vxy ;-v_uxy])

% Get u(x,y) as a matrix of coefficients
uxy = GetAsMatrix(v_uxy, m1-t1, m2-t2);

% Get v(x,y) as a matrix of coefficients
vxy = GetAsMatrix(v_vxy,n1-t1,n2-t2);


switch SETTINGS.SYLVESTER_BUILD_METHOD
    case 'T'
        
        vxy = GetWithoutBinomials_Bivariate(vxy);
        uxy = GetWithoutBinomials_Bivariate(uxy);

    case 'DT'
        
        vxy = GetWithoutBinomails_Bivariate(vxy);
        uxy = GetWithoutBinomials_Bivariate(uxy);
        
    case 'DTQ'
        
    case 'TQ'
        
        

    otherwise
        error('err')
end

end