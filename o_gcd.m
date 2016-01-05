function [dxy_matrix_calc] = o_gcd(ex_num,el,BOOL_PREPROC,BOOL_SNTLN,SEED)
%%
% Given an example number and set of parameters, obtain GCD of the two
% polynomials fxy and gxy in the given example file. 
% Where the polynomials f(x,y) and g(x,y) are defined as polynomails in 
% the Bernstein Basis.

% %                             Inputs 

% ex_num - Example Number

% el - Lower noise level

% BOOL_PREPROC ('y'/'n')
%       y - Include Preprocessing
%       n - Exclude Preprocessing

% BOOL_SNTLN ('y'/'n')
%       y - Include SNTLN
%       n - Exclude SNTLN

% SEED - SEED Number for noise generation


%%
%                           Set Variables

% bool_Q ('y'/'n')
%       y - Include the matrix Q in the Sylvester matrix such that it is
%       given by D^{-1}T(f,g)Q
%       n - Exclude the matrix Q from the Sylvester matrix such that it is
%       given by D^{-1}T(f,g)
global bool_Q
bool_Q = input('Do you wish to include Q in the coefficient matrix (y)es/(n)o :  ','s');

% bool_preproc ('y'/'n')
%       y - Include the three preprocessing operations.
%       n - Exclude the three preprocessing operations.
global bool_preproc
bool_preproc = BOOL_PREPROC;


% bool_noise ('y'/'n')
%       y - Add noise to the coefficients of input polynomial f
%       n - Use exact form of input polynomial f
global bool_noise
bool_noise = 'y';

% bool_SNTLN ('y'/'n')
%       y - Include SNTLN
%       n - Exclude SNTLN
global bool_SNTLN
bool_SNTLN = BOOL_SNTLN;

% bool_plotgraphs ('y'/'n')
%       y - plot 
%       n - 
global bool_plotgraphs
bool_plotgraphs = 'n';

% seed - SEED Number for noise generation
global seed
seed = SEED;

global threshold
threshold = 1;
%%
%                   Get Example


[fxy_matrix_exact, gxy_matrix_exact,...
    uxy_matrix_exact,vxy_matrix_exact,...
    dxy_matrix_exact,...
    m,m1,m2,...
    n,n1,n2,...
    t_exact,t1_exact,t2_exact] = Examples_GCD(ex_num);




fprintf('\n')
fprintf('----------------------------------------------------------------\n')
fprintf('Input Polynomials Degrees:\n')
fprintf('m  : %i \n',m)
fprintf('m1 : %i \n',m1)
fprintf('m2 : %i \n\n',m2)
fprintf('n  : %i \n',n)
fprintf('n1 : %i \n',n1)
fprintf('n2 : %i \n\n',n2)

fprintf('t  : %i \n',t_exact)
fprintf('t1 : %i \n',t1_exact)
fprintf('t2 : %i \n',t2_exact)
fprintf('----------------------------------------------------------------\n')
fprintf('\n')


% Add Noise to the coefficients
switch bool_noise
    case 'y'
        % Add noise to the coefficients of f and g
        [fxy_matrix, noise_mat_f] = Noise2(fxy_matrix_exact,el);
        [gxy_matrix, noise_mat_g] = Noise2(gxy_matrix_exact,el);
    case 'n'
    otherwise 
        error('noise value either y or n')
end


% Plot the surfaces of the two polynomials fxy and gxy
if bool_plotgraphs == 'y'
    plot_fxy_gxy(fxy_matrix,gxy_matrix);
end

%%

%                   Calculate GCD

% Calculate the gcd, and quotient polynomials of f(x,y) and g(x,y)
[uxy_matrix_calc, vxy_matrix_calc, dxy_matrix_calc] = o1(fxy_matrix,gxy_matrix,...
    m,n);

%% 

%                   Results.


fprintf('----------------------------------------------------------------')
fprintf('\n')
fprintf('Compare Exact Coefficients with Computed Coefficients of u(x,y):')
uxy_matrix_calc ./ uxy_matrix_calc(1,1)
uxy_matrix_exact./uxy_matrix_exact(1,1)

(uxy_matrix_calc ./ uxy_matrix_calc(1,1)) ./ (uxy_matrix_exact./uxy_matrix_exact(1,1))

fprintf('----------------------------------------------------------------')
fprintf('\n')
fprintf('Compare Exact Coefficients with Computed Coefficients of v(x,y): ')
vxy_matrix_calc ./ vxy_matrix_calc(1,1)
vxy_matrix_exact./vxy_matrix_exact(1,1)
(vxy_matrix_calc ./ vxy_matrix_calc(1,1)) ./ (vxy_matrix_exact./vxy_matrix_exact(1,1))

fprintf('----------------------------------------------------------------')
fprintf('\n')

fprintf('Compare Exact coefficients with computed coefficients of d(x,y)')
dxy_matrix_calc./dxy_matrix_calc(1,1)
dxy_matrix_exact./dxy_matrix_exact(1,1)
(dxy_matrix_exact./dxy_matrix_exact(1,1)) ./ (dxy_matrix_calc./dxy_matrix_calc(1,1))

% Get Euclidean Distance between the exact and calculated d(x,y)
EuclideanDistance_d = ...
norm(abs(...
    (dxy_matrix_exact./dxy_matrix_exact(1,1)) ...
    - ...
    dxy_matrix_calc./dxy_matrix_calc(1,1)));

fprintf('Euclidean distance of d(x,y) exact to d(x,y) computed is \n %0.6e \n',EuclideanDistance_d)
fprintf('\n')
fprintf('----------------------------------------------------------------')
fprintf('\n')
end