function [uxy_matrix_calc, vxy_matrix_calc] = GetCofactors(fxy_matrix, gxy_matrix, t1,t2)
% Given two input polynomials and the degree of the GCD return the quotient
% polynomials u(x,y) and v(x,y)
%
%   Inputs.
%
%   fxy_matrix : Coefficients of polynomial f(x,y)
%
%   gxy_matrix : Coefficients of polynomial g(x,y)
%
%   t1 : Degree of d(x,y) with respect to x
%
%   t2 : Degree of d(x,y) with respect to y
%

% Initialise Global Variables
global SETTINGS

%%
% Get the degrees of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get the degrees of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% Build the (t1,t2)-th subresultant


Sk1k2 = BuildDTQ(fxy_matrix, gxy_matrix,t1,t2);


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

num_coeff_v = (n1-t1+1) * (n2-t2+1);
num_coeff_u = (m1-t1+1) * (m2-t2+1);


% get coefficients of u(w,w) and v(w,w)
vww_calc = vecx(1:num_coeff_v);
uww_calc = -vecx(num_coeff_v+1:end);


%% Get the value of fv-gu
norm(Sk1k2 * [vww_calc ;-uww_calc])

%% Obtain u(x,y) in its matrix form
% Arrange uw into a matrix form based on its dimensions.
uxy_matrix_calc = GetAsMatrix(uww_calc,m1-t1,m2-t2);

%% Obtain v(x,y) in its matrix form
% Arrange vw into a matrix form based on their dimensions.
vxy_matrix_calc = GetAsMatrix(vww_calc,n1-t1,n2-t2);


% If we excluded Q from the coefficient matrix, then remove the binomial 
% coefficients from v(x,y) and u(x,y)

switch SETTINGS.BOOL_Q
    case 'y'
        % Do nothing
    case 'n'
        % If Q is not included in the Sylvester matrix, then Q is included
        % in the solution vector, so strip binomial coefficients from
        % v(x,y) and u(x,y)
        vxy_matrix_calc = GetWithoutBinomails(vxy_matrix_calc);
        
        uxy_matrix_calc = GetWithoutBinomials(uxy_matrix_calc);

    otherwise
        error('err')
end

end