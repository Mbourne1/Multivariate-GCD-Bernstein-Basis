function [dxy_matrix_calc] = o_gcd(ex_num,el,bool_preproc,low_rank_approx_method)
% o_gcd(ex_num,el,bool_preproc,bool_sntln)
%
% Given an example number and set of parameters, obtain GCD of the two
% polynomials f(x,y) and g(x,y) in the given example file. 
% Where the polynomials f(x,y) and g(x,y) are defined as polynomails in 
% the Bernstein Basis.
%
% %                             Inputs 
%
% ex_num - Example Number
%
% el - Lower noise level
%
% bool_preproc ('y'/'n')
%       'y' : Include Preprocessing
%       'n' : Exclude Preprocessing
%
% low_rank_approx_method :
%       'Standard SNTLN' : Include SNTLN
%       'None' : Exclude SNTLN
%

% %
% Set Variables
SetGlobalVariables()


% %
% Get Example

[fxy_matrix_exact, gxy_matrix_exact,...
    uxy_matrix_exact,vxy_matrix_exact,...
    dxy_matrix_exact,...
    m,m1,m2,...
    n,n1,n2,...
    t_exact,t1_exact,t2_exact] = Examples_GCD(ex_num);

display(fxy_matrix_exact)
display(gxy_matrix_exact)
display(dxy_matrix_exact)



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

% %
% Add Noise to the coefficients
switch BOOL_NOISE
    case 'y'
        % Add noise to the coefficients of f and g
        [fxy_matrix, ~] = Noise2(fxy_matrix_exact,el);
        [gxy_matrix, ~] = Noise2(gxy_matrix_exact,el);
    case 'n'
    otherwise 
        error('noise value either y or n')
end


% Plot the surfaces of the two polynomials fxy and gxy
if PLOT_GRAPHS == 'y'
    plot_fxy_gxy(fxy_matrix,gxy_matrix);
end

% % Calculate GCD

% Calculate the gcd, and quotient polynomials of f(x,y) and g(x,y)
[uxy_matrix_calc, vxy_matrix_calc, dxy_matrix_calc] = o1(fxy_matrix,gxy_matrix,...
    m,n);

% % Results.

PrintoutCoefficients('u',uxy_matrix_calc,uxy_matrix_exact)
PrintoutCoefficients('v',vxy_matrix_calc,vxy_matrix_exact)
PrintoutCoefficients('d',dxy_matrix_calc,dxy_matrix_exact)

end

function [] = PrintoutCoefficients(u,uxy_matrix_calc,uxy_matrix_exact)
fprintf('----------------------------------------------------------------')
fprintf('\n')
fprintf('Compare Exact Coefficients with Computed Coefficients of %s(x,y):',u)
uxy_matrix_calc ./ uxy_matrix_calc(1,1)
uxy_matrix_exact./uxy_matrix_exact(1,1)

end
