
function Sk = BuildSubresultant(fxy_matrix,gxy_matrix,k1,k2,th1,th2)
% Build the sylvester subresultant matrix S_{k1,k2}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                         Inputs.


% fxy_matrix_bi - coefficients of polynomial fxy in the scaled bernstein
% basis (including the binomial coefficients) a_{i,j}
% \binom{m_{1}}{i}\binom{m_{2}}{j}

% gxy_matrix_bi - coefficients of polynomial gxy in the scaled bernstein
% basis (including the binomial coefficients) b_{i,j}
% \binom{n_{1}}{i}\binom{n_{2}}{j}

% k1 - the degree k_{1} with respect to x of the polynomial d_{k_{1},k_{2}}

% k2 - The degree k_{2} with respect to y of the polynomial d_{k_{1},k_{2}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                     Global Variables


% bool_Q : Global boolean value Q,
%   1 - Include Q in the Sylvester Subresultant matrix D^{-1}T(f,g)Q
%   0 - Exclude Q from the Sylvester Subresultant matrix D^{-1}T(f,g)

global bool_Q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %                         Code

[rows,cols] = size(fxy_matrix);

% get the degrees of fxy with respect to x and y
m1 = rows - 1;
m2 = cols - 1;

% get the degrees of gxy with respect to x and y
[rows,cols] = size(gxy_matrix);
n1 = rows - 1;
n2 = cols - 1;

% Build two Cauchy matrices, the first for coefficients of fxy and the
% second for the coefficients of gxy
C_f = BuildT1(fxy_matrix,n1,n2,k1,k2,th1,th2);
C_g = BuildT1(gxy_matrix,m1,m2,k1,k2,th1,th2);

% Build the diagonal matrix D^{-1}
D = BuildD(k1,k2,m1,m2,n1,n2);



switch bool_Q
    case 1
        % Build the diagonal matrix Q such that Q * [v \\ u] gives the
        % coefficients of u and v in the scaled bernstein basis
        Q = BuildQ(k1,k2,m1,m2,n1,n2);
        Sk = D*[C_f,C_g]*Q;
    case 0
        Sk = D*[C_f,C_g];
end
end
