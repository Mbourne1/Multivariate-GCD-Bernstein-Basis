function DYQ = BuildDYG_SNTLN(m1,m2,n1,n2,k1,k2,x,alpha,th1,th2)
% BuildY(m1,m2,n1,n2,k1,k2,idx_col,x_ls,alpha,th1,th2)
%
%
% Build the matrix Y such that D*Y(x)*G*[f;g] = A_{t}(f,g)*x, 
% where A_{t}(f,g) = (D*T(f,g)*Q)*M where M removes an optimal column of
% T(f,g).
% This function is used in STLN()
%
% Inputs.
%
% m1 : Degree of polynomial f(x,y) with respect to x
%
% m2 : Degree of polynomial f(x,y) with respect to y
%
% n1 : Degree of polynomial g(x,y) with respect to x
%
% n2 : Degree of polynomial g(x,y) with respect to y
%
% k1 : Degree of polynomial d(x,y) with respect to x
%
% k2 : Degree of polynomial d(x,y) with respect to y
%
% x : 
%
% alpha :
%
% th1 :
%
% th2 :
%
% % Outputs
%
% DYQ :

% Get degree of polynomial u(x,y)
m1_k1 = m1 - k1;
m2_k2 = m2 - k2;

% Get degree of polynomial v(x,y)
n1_k1 = n1 - k1; 
n2_k2 = n2 - k2;


% Separate the x into coefficients of u and coefficients of v
% The first (n1-t1+1) x (n2-t2+1) coefficients are of v
nCoeffs_x1 = (n1_k1+1) * (n2_k2+1);

% Get the vector of coefficients of v(w,w)
x1_vec = x(1:nCoeffs_x1);

% Get the vector of coefficients of u(w,w)
x2_vec = x(nCoeffs_x1+1:end);

% Get the vector x1(w,w) and x2(w,w) as a matrix
x1_ww = GetAsMatrix(x1_vec,n1_k1,n2_k2);
x2_ww = GetAsMatrix(x2_vec,m1_k1,m2_k2);

% Build the matrix D^{-1}_{m1+n1-k1,m2+n2-k2}
D = BuildD(m1,m2,n1-k1,n2-k2);

% Construct T1(x1) and T2(x2).
Y1 = BuildT1(x1_ww,m1,m2);
Y2 = BuildT1(x2_ww,n1,n2);

% Build the matrices Q_{m1,m2} and Q_{n1,n2}
Q1 = BuildQ1(m1,m2);
Q2 = BuildQ1(n1,n2);

% Get the matrix of thetas corresponding to (and exlcluding) coefficients
% of f(x,y)
th_fww = GetAsVector(GetWithThetas(ones(m1+1,m2+1),th1,th2));
th_fww = diag(th_fww);

% Get the matrix of thetas corresponding to (and exlcluding) coefficients
% of g(x,y)
th_gww = GetAsVector(GetWithThetas(ones(n1+1,n2+1),th1,th2));
th_gww = diag(th_gww);


% multiply by the alpha of g
DYQ = [D*Y1*Q1*th_fww alpha.*D*Y2*Q2*th_gww];
end