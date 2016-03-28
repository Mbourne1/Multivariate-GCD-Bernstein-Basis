function Sk = BuildDTQ(fxy_matrix,gxy_matrix,k1,k2)
% BuildSubresultant(fxy_matrix_n,gxy_matrix_n,k1,k2,alpha,th1,th2)
%
% Build the sylvester subresultant matrix S_{k1,k2}.
%
% Inputs.
%
%
% fxy_matrix_n : Coefficients of the polynomial f(x,y)
%
% gxy_matrix_n : Matrix of coefficients of polynomial gxy in the scaled 
%                bernstein basis (including the binomial coefficients) b_{i,j}
%                \binom{n_{1}}{i}\binom{n_{2}}{j}
%
% k1 : The degree k_{1} with respect to x of the polynomial d_{k_{1},k_{2}}
%
% k2 : The degree k_{2} with respect to y of the polynomial d_{k_{1},k_{2}}
%
%

% %                     Global Variables


% BOOL_Q : Global boolean value Q,
%   1 - Include Q in the Sylvester Subresultant matrix D^{-1}T(f,g)Q
%   0 - Exclude Q from the Sylvester Subresultant matrix D^{-1}T(f,g)

global BOOL_Q

%

% %                         Code

% Get the degree of f(x,y) with respect to x and y
[m1,m2] = GetDegree(fxy_matrix);

% Get the degree of g(x,y) with respect to x and y
[n1,n2] = GetDegree(gxy_matrix);

% Build two Cauchy matrices, the first for coefficients of fxy and the
% second for the coefficients of gxy
C_f = BuildT1(fxy_matrix,n1-k1,n2-k2);
C_g = BuildT1(gxy_matrix,m1-k1,m2-k2);

% Build the diagonal matrix D^{-1}
D = BuildD(k1,k2,m1,m2,n1,n2);



switch BOOL_Q
    case 'y'
        % Build the diagonal matrix Q such that Q * [v \\ u] gives the
        % coefficients of u and v in the scaled bernstein basis
        Q = BuildQ(n1-k1,n2-k2,m1-k1,m2-k2);
        Sk = D*[C_f, C_g]*Q;
    case 'n'
        
        Sk = D*[C_f, alpha.* C_g];
end



end
