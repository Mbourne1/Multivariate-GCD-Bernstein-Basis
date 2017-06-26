function DPG = BuildDPG_SNTLN(m1,m2,n1,n2,th1,th2,alpha,k1,k2,idx_col)
% BuildP(m1,m2,n1,n2,th1,th2,alpha,k1,k2,idx_col)
%
% Build the matrix P, which is used in the SNTLN function, where a column
% of the Sylvester matrix S_{k_{1},k_{2}} whose index is given, can be
% described as P*[f;g].
%
% Inputs
%
% m1 : (Int) Degree of polynomial f(x,y) with respect to x
%
% m2 : (Int) Degree of polynomial f(x,y) with respect to y
%
% n1 : (Int) Degree of polynomial g(x,y) with respect to x
%
% n2 : (Int) Degree of polynomial g(x,y) with respect to y
%
% th1 : (Float) Optimal value of \theta_{1}
%
% th2 : (Float) Optimal value of \theta_{2}
%
% k1 : (Int) Degree of polynomial d(x,y) with respect to x
%
% k2 : (Int) Degree of polynomial d(x,y) with respect to y
%
% idx_col : (Int) Index of column removed from the Sylvester Subresultant
% matrix S_{k_{1},k_{2}}(f,g) such that A_{k_{1},k_{2}}x = c_{k_{1},k_{2}}.
%
% Outputs.
%
% P : (Matrix) P, such that P*[f;g] = c_{k_{1},k_{2}}


% Get the number of columns in the first partition of the Sylvester
% subresultant S_{t_{1},t_{2}}.
nColumns_T1 = (n1 - k1 + 1) * (n2 - k2 + 1);

% Get the number of coefficients in the polynomial f(x,y)
nCoefficients_fxy = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
nCoefficients_gxy = (n1+1) * (n2+1);

% Get the number of rows in the Sylvester matrix S_{k_{1},k_{2}}
nRows_Sk1k2 = (m1+n1-k1+1)*(m2+n2-k2+1);

% If column is in the first partition of S_{k_{1},k_{2}}
if idx_col <= nColumns_T1
    
    % Build the matrix P1
    P1 = BuildP1(m1, m2, n1, n2, k1, k2, idx_col);
    
    % Build the matrix P2
    P2 = zeros(nRows_Sk1k2, nCoefficients_gxy);
    
    
else % Column is in the second partiton of S_{k_{1},k_{2}}
    
    
    % Build the matrix P1
    P1 = zeros(nRows_Sk1k2, nCoefficients_fxy);
    
    % Get the index of the column with respect to T(g) only.
    idx_col_rel = idx_col - nColumns_T1;
    P2 = BuildP1(n1, n2, m1, m2, k1, k2, idx_col_rel);
    
    
end

% Build the matrix D_{m1+n1-k1,m2+n2-k2}^{-1}
D = BuildD_Bivariate_2Polys(m1, m2, n1-k1, n2-k2);

% Build the matrix Q_{m}
Q_m = BuildQ1_Bivariate(m1, m2);

%Build the matrix Q_{n}
Q_n = BuildQ1_Bivariate(n1, n2);

% Get the thetas corresponding to coefficients of f(x,y)
th_f = GetAsVector_Version1(GetWithThetas(ones(m1 + 1, m2 + 1), th1, th2));
th_f = diag(th_f);

% Get the thetas corresponding to coefficients of g(x,y)
th_g = GetAsVector_Version1(GetWithThetas(ones(n1 + 1, n2 + 1), th1, th2));
th_g = diag(th_g);

% Construct the matrix P
DPG =  D*[P1*th_f*Q_m alpha.*P2*th_g*Q_n];

end