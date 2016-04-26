function [dxy_matrix_calc] = o_gcd(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
% o_gcd(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
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
% mean_method
%       'None'
%       'Geometric Mean Matlab Method'
%
% bool_preproc ('y'/'n')
%       'y' : Include Preprocessing
%       'n' : Exclude Preprocessing
%
% low_rank_approx_method :
%       'Standard SNTLN' : Include SNTLN
%       'Standard STLN : Include STLN
%       'None' : Exclude SNTLN
%

% %
% Set Variables
SetGlobalVariables(mean_method,bool_alpha_theta,low_rank_approx_method)


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
global BOOL_NOISE
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
global PLOT_GRAPHS
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

dxy_error = GetDistance('d',dxy_matrix_calc,dxy_matrix_exact);

PrintToFile(m,n,dxy_error);


end


function [] = PrintoutCoefficients(name,matrix_calc,matrix_exact)
fprintf('----------------------------------------------------------------')
fprintf('\n')
fprintf('Compare Exact Coefficients with Computed Coefficients of %s(x,y):',name)

matrix_calc  = normalise(matrix_calc);
matrix_exact = normalise(matrix_exact);

display(matrix_calc)
display(matrix_exact)


end

function [dist] = GetDistance(name,matrix_calc,matrix_exact)

matrix_calc = normalise(matrix_calc);
matrix_exact = normalise(matrix_exact);

fprintf('Analysis of Coefficients of %s(x,y) computed vs %s(x,y) exact: \n',name,name)
fprintf('Distance between exact and calculated matrix:')
dist = (norm(matrix_exact,'fro') - norm(matrix_calc,'fro') )./ norm(matrix_exact,'fro');


end


function []= PrintToFile(m,n,error_dx)

global NOISE
global BOOL_ALPHA_THETA

fullFileName = 'o_gcd_results.txt';


if exist('o_gcd_results.txt', 'file')
    fileID = fopen('o_gcd_results.txt','a');
    fprintf(fileID,'%5d \t %5d \t %5d \t %s \t %s \t %s\n',...
        m,n,error_dx,BOOL_ALPHA_THETA, NOISE);
    fclose(fileID);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
  uiwait(msgbox(warningMessage));
end

end
