function [dxy_calc] = o_gcd(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method)
% o_gcd(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
%
% Given an example number and set of parameters, obtain GCD of the two
% polynomials f(x,y) and g(x,y) in the given example file.
% Where the polynomials f(x,y) and g(x,y) are defined as polynomails in
% the Bernstein Basis.
%
% % Inputs
%
% ex_num - Example Number
%
% el - Lower noise level
%
% em - Upper noise level
%
% mean_method
%       'None'
%       'Geometric Mean Matlab Method'
%
% bool_alpha_theta ('y'/'n')
%       'y' : Include Preprocessing
%       'n' : Exclude Preprocessing
%
% low_rank_approx_method :
%       'Standard SNTLN' : Include SNTLN
%       'Standard STLN : Include STLN
%       'None' : Exclude SNTLN
%
% % Examples
%
% >> o_gcd('1',1e-12,1e-10,'Geometric Mean Matlab Method','y','None')

% %
% Set Variables
SetGlobalVariables(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method)
global SETTINGS

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

DegreeStructure()


% %
% Add Noise to the coefficients


% Add noise to the coefficients of f and g
[fxy_matrix, ~] = Noise2(fxy_matrix_exact,emin,emax);
[gxy_matrix, ~] = Noise2(gxy_matrix_exact,emin,emax);



% Plot the surfaces of the two polynomials fxy and gxy

if SETTINGS.PLOT_GRAPHS == 'y'
    plot_fxy_gxy(fxy_matrix,gxy_matrix);
end

% % Calculate GCD
lower_limit = 1;
upper_limit = min(m,n);

% Calculate the gcd, and quotient polynomials of f(x,y) and g(x,y)
[fxy_calc,gxy_calc,dxy_calc,uxy_calc, vxy_calc,t,t1,t2] = o_gcd_mymethod(fxy_matrix,gxy_matrix,...
    m,n,[lower_limit,upper_limit]);

% % Results.

PrintoutCoefficients('u',uxy_calc,uxy_matrix_exact)
PrintoutCoefficients('v',vxy_calc,vxy_matrix_exact)
PrintoutCoefficients('d',dxy_calc,dxy_matrix_exact)

error.dxy = GetDistance('d',dxy_calc,dxy_matrix_exact);
error.uxy = GetDistance('u',uxy_calc,uxy_matrix_exact);
error.vxy = GetDistance('v',vxy_calc,vxy_matrix_exact);
PrintToFile(m,n,error);


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

fprintf([mfilename sprintf('Analysis of Coefficients of %s(x,y) computed vs %s(x,y) exact: \n',name,name)])
fprintf([mfilename sprintf('Distance between exact and calculated matrix: \n')])
dist = (norm(matrix_exact,'fro') - norm(matrix_calc,'fro') )./ norm(matrix_exact,'fro');


end


function []= PrintToFile(m,n,error)

global SETTINGS


fullFileName = 'Results_o_gcd.txt';


if exist('Results_o_gcd.txt', 'file')
    fileID = fopen('Results_o_gcd.txt','a');
    fprintf(fileID,'%s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s \n',...
        SETTINGS.EX_NUM,...
        num2str(m),...
        num2str(n),...
        error.dxy,...
        error.uxy,...
        error.vxy,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.EMIN);
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end

end
