function [uxy, vxy] = GetCofactors_Bivariate_2Polys(fxy, gxy, t1, t2)
% Given two input polynomials and the degree of the GCD return the quotient
% polynomials u(x,y) and v(x,y)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% t1 : (Int) Degree of d(x,y) with respect to x
%
% t2 : (Int) Degree of d(x,y) with respect to y
%
% % Outputs
%
% uxy : (Matrix) Coefficients of polynomial u(x,y), the quotient polynomial
% such that f(x,y)/u(x,y) = d(x,y)
%
% vxy : (Matrix) Coefficients of polynomial v(x,y), the quotient polynomial
% such that g(x,y)/v(x,y) = d(x,y)


% Initialise Global Variables
global SETTINGS

% Get the degrees of polynomial f(x,y) and g(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Build the (t1,t2)-th subresultant matrix
Sk1k2 = BuildSubresultant_Bivariate_2Polys(fxy, gxy, t1, t2);

% Find Optimal column for removal from St
% given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized
idxOptimalColumn = GetOptimalColumn(Sk1k2);

% Get the matrix A_{t_{1},t_{2}}
At = Sk1k2;
At(:, idxOptimalColumn) = [];

% Get the vector c_{t_{1},t_{2}} removed from S_{t_{1},t_{2}}
ct = Sk1k2(:, idxOptimalColumn);


% Get the coefficients for u(x,y) and v(x,y)
x_ls = SolveAx_b(At,ct);


% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1: (idxOptimalColumn) -1);
    -1;
    x_ls(idxOptimalColumn : end);
    ];

% Get number of coefficients in v(x,y) and u(x,y)
nCoefficients_vxy = (n1 - t1 + 1) * (n2 - t2 + 1);
nCoefficients_uxy = (m1 - t1 + 1) * (m2 - t2 + 1);

% Get coefficients of u(x,y) as a vector
v_vxy = vecx(1:nCoefficients_vxy);

% Get coefficients of v(w,w) as a vector
v_uxy = -1 .* vecx(nCoefficients_vxy+1 : nCoefficients_uxy + nCoefficients_vxy);

% Get u(x,y) and v(x,y) as matrices of coefficients
uxy = GetAsMatrix(v_uxy, m1 - t1, m2 - t2);
vxy = GetAsMatrix(v_vxy, n1 - t1, n2 - t2);


switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    case 'T'
        
        % Get u(x,y) and v(x,y) as a matrix of coefficients
        uxy = GetAsMatrix(v_uxy, m1 - t1, m2 - t2);
        vxy = GetAsMatrix(v_vxy, n1 - t1, n2 - t2);
        vxy = GetWithoutBinomials_Bivariate(vxy);
        uxy = GetWithoutBinomials_Bivariate(uxy);
        
    case 'DT'

        % Get u(x,y) and v(x,y) as a matrix of coefficients
        uxy = GetAsMatrix(v_uxy, m1 - t1, m2 - t2);
        vxy = GetAsMatrix(v_vxy, n1 - t1, n2 - t2);
        vxy = GetWithoutBinomials_Bivariate(vxy);
        uxy = GetWithoutBinomials_Bivariate(uxy);
        
    case 'DTQ'

        % Get u(x,y) and v(x,y) as a matrix of coefficients
        uxy = GetAsMatrix(v_uxy, m1 - t1, m2 - t2);
        vxy = GetAsMatrix(v_vxy, n1 - t1, n2 - t2);
        
    case 'DTQ Version 2'
        
        % Get u(x,y) as a matrix of coefficients
        uxy = GetAsMatrix_Version2(v_uxy, m1 - t1, m2 - t2);
        
        % Get v(x,y) as a matrix of coefficients
        vxy = GetAsMatrix_Version2(v_vxy, n1 - t1, n2 - t2);
        
        
    case 'TQ'
        
        % Get u(x,y) and v(x,y) as a matrix of coefficients
        uxy = GetAsMatrix(v_uxy, m1 - t1, m2 - t2);
        vxy = GetAsMatrix(v_vxy, n1 - t1, n2 - t2);
        
    case 'DTQ Denominator Removed'
        
        % Get u(x,y) and v(x,y) as a matrix of coefficients
        uxy = GetAsMatrix(v_uxy, m1 - t1, m2 - t2);
        vxy = GetAsMatrix(v_vxy, n1 - t1, n2 - t2);
        
    otherwise
        error('err')
end

end