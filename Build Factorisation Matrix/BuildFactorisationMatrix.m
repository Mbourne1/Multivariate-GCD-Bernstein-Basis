function C_uv = BuildFactorisationMatrix(uxy, vxy, t1, t2)
% Build the matrix C consisting of two convolution matrices C(u) and C(v)
% stacked on top of one another C = [ C(u) ; C(v)]. This matrix is used in
% the computation of the coefficients of the GCD d(x,y) 
%
%
% % Inputs
%
% uxy : (Matrix) Coefficients of polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of polynomial v(x,y)
%
% t1 : (Int) Degree of common divisor with respect to x
%
% t2 : (Int) Degree of common divisor with respect to y
%
% % Outputs
%
% C_uv : (Matrix) 

global SETTINGS

% Get degree of u(x,y) and v(x,y)
[m1_t1, m2_t2] = GetDegree_Bivariate(uxy);
[n1_t1, n2_t2] = GetDegree_Bivariate(vxy);

% Get degree structure of f(x,y)
m1 = m1_t1 + t1;
m2 = m2_t2 + t2;

% Get degree structure of g(x,y)
n1 = n1_t1 + t1;
n2 = n2_t2 + t2;

% Build matrix H
H = BuildH_Bivariate(m1, m2, n1, n2);

% Build the matrix C_{}(u(x,y))
C1_u = BuildT1_Bivariate(uxy, t1, t2);

% Build Matrix C_{}(v(x,y)
C2_v = BuildT1_Bivariate(vxy, t1, t2);

C = [ ...
    C1_u;
    C2_v ];

% Buid matrix G
G = BuildG_Bivariate(t1, t2);

switch SETTINGS.FACTORISATION_BUILD_METHOD

    case 'HCG'

        % Build the Coefficient Matrix HCG 
        C_uv = H*C*G;

    case 'CG'
        
        C_uv = C*G;
        
    case 'HC'

        C_uv = H*C;
        
    otherwise
        
        error('factorisation method not valid')
        
end