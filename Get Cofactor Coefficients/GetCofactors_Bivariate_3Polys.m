function [uxy, vxy, wxy] = GetCofactors_Bivariate_3Polys(fxy, gxy, hxy, t1, t2)
% Given two input polynomials and the degree of the GCD return the quotient
% polynomials u(x,y) and v(x,y)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% wxy : (Matrix) Coefficients of polynomial h(x,y)
%
% t1 : (Int) Degree of d(x,y) with respect to x
%
% t2 : (Int) Degree of d(x,y) with respect to y
%
% % Outputs
%
% uxy : (Matrix) Coefficients of polynomial u(x,y), the quotient polynomial such
% that f(x,y)/u(x,y) = d(x,y)
%
% vxy : (Matrix) Coefficients of polynomial u(x,y), the quotient polynomial such
% that g(x,y)/v(x,y) = d(x,y)
%
% wxy : (Matrix) Coefficients of polynomial u(x,y), the quotient polynomial such
% that h(x,y)/w(x,y) = d(x,y)



% Initialise Global Variables
global SETTINGS


% Get the degrees of polynomial f(x,y), g(x,y) and h(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);

% Build the (t1,t2)-th subresultant
Sk1k2 = BuildSubresultant_Bivariate_3Polys(fxy, gxy, hxy, t1, t2);


% Find Optimal column for removal from S_{t1,t2}
% Given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized
opt_col_index = GetOptimalColumn(Sk1k2);

% Get the matrix A_{t_{1},t_{2}} 
At = Sk1k2;
At(:, opt_col_index) = [];

% Get the vector c_{t_{1},t_{2}} removed from S_{t_{1},t_{2}}
ct = Sk1k2(:, opt_col_index);


% % Get the coefficients for u(x,y) and v(x,y)
x_ls = SolveAx_b(At,ct);


% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col_index)-1);
    -1;
    x_ls(opt_col_index:end);
    ];

nCoefficients_vxy = (n1 - t1 + 1) * (n2 - t2 + 1);

nCoefficients_wxy = (o1 - t1 + 1) * (o2 - t2 + 1);


% Get coefficients of u(x,y) as a vector 
v_vxy = vecx(1:nCoefficients_vxy);

v_wxy = vecx(nCoefficients_vxy + 1 : nCoefficients_vxy + nCoefficients_wxy);

% Get coefficients of v(w,w) as a vector
v_uxy = -1 .* vecx(nCoefficients_vxy + nCoefficients_wxy + 1 : end);


% Get u(x,y), v(x,y) and w(x,y) as matrices of coefficients
uxy = GetAsMatrix(v_uxy, m1 - t1, m2 - t2);
vxy = GetAsMatrix(v_vxy, n1 - t1, n2 - t2);
wxy = GetAsMatrix(v_wxy, o1 - t1, o2 - t2);

switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    case 'T'
              
        vxy = GetWithoutBinomials_Bivariate(vxy);
        uxy = GetWithoutBinomials_Bivariate(uxy);
        wxy = GetWithoutBinomials_Bivariate(wxy);
        
    case 'DT'
        
        vxy = GetWithoutBinomials_Bivariate(vxy);
        uxy = GetWithoutBinomials_Bivariate(uxy);
        wxy = GetWithoutBinomials_Bivariate(wxy);
        
    case 'DTQ'
        
    case 'TQ'
        
    case 'DTQ Denominator Removed'

    otherwise
        error('err %s is not a valid SYLVESTER BUILD METHOD', SETTINGS.SYLVESTER_MATRIX_VARIANT)
end

end