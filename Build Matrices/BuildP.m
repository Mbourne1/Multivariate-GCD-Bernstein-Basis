function P = BuildP(m1,m2,n1,n2,th1,th2,alpha,k1,k2,idx_col)
% BuildP(m1,m2,n1,n2,th1,th2,alpha,k1,k2,idx_col)
%
% Build the matrix P, which is used in the SNTLN function, where a column
% of the Sylvester matrix S_{k_{1},k_{2}} whose index is given, can be
% described as P*[f;g].
%
% Inputs
%
% m1 : Degree of polynomial f(x,y) with respect to x
%
% m2 : Degree of polynomial f(x,y) with respect to y
%
% n1 : Degree of polynomial g(x,y) with respect to x
%
% n2 : Degree of polynomial g(x,y) with respect to y
%
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}
%
% k1 : Degree of polynomial d(x,y) with respect to x
%
% k2 : Degree of polynomial d(x,y) with respect to y
%
% idx_col : Index of column removed from the Sylvester Subresultant
% matrix S_{k_{1},k_{2}}(f,g) such that A_{k_{1},k_{2}}x = c_{k_{1},k_{2}}.
%
% Outputs.
%
% P : Matrix P, such that P*[f;g] = c_{k_{1},k_{2}}


% Get the number of columns in the first partition of the Sylvester
% subresultant S_{t_{1},t_{2}}.
num_cols_T1 = (n1 - k1 + 1) * (n2 - k2 + 1);

% Get the number of coefficients in the polynomial f(x,y)
num_coeff_f = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
num_coeff_g = (n1+1) * (n2+1);

% Get the number of rows in the Sylvester matrix S_{k_{1},k_{2}}
nRows = (m1+n1-k1+1)*(m2+n2-k2+1);

% If column is in the first partition of S_{k_{1},k_{2}}
if idx_col <= num_cols_T1
    
    % Build the matrix P1
    P1 = BuildP1(m1,m2,n1,n2,th1,th2,k1,k2,idx_col);
    
    % Build the matrix P2
    P2 = zeros(nRows,num_coeff_g);
    
    
else % Column is in the second partiton of S_{k_{1},k_{2}}
    
    
    % Build the matrix P1
    P1 = zeros(nRows,num_coeff_f);
    
    % Get the index of the column with respect to T(g) only.
    idx_col_rel = idx_col - num_cols_T1;
    P2 = BuildP1(n1,n2,m1,m2,th1,th2,k1,k2,idx_col_rel);
    
    
end

% Construct the matrix P
P =  [P1 alpha.*P2];

end