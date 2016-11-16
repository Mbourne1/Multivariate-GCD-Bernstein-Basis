function [dxy_calc] = o_gcd(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method)
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
% low_rank_approx_method 
%       'Standard SNTLN' : Include SNTLN
%       'Standard STLN : Include STLN
%       'None' : Exclude SNTLN
%
% apf_method
%       'None'
%       'Standard APF Nonlinear'
%       'Standard APF Linear'
%
% sylvester_build_method 
%       'T'
%       'DT'
%       'DTQ'
%       'TQ'
%
% % Examples
%
% >> o_gcd('1',1e-12,1e-10,'Geometric Mean Matlab Method','y','Standard STLN','None','DTQ')

% %
% Set Variables

SetGlobalVariables(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method,sylvester_build_method,apf_method)

% Add subfolders
restoredefaultpath

addpath(...
    'Basis Conversion',...
    'Bernstein Methods',...
    'Build Matrices',...
    'Formatting',...
    'Get Cofactors',...
    'Get GCD Coefficients',...
    'Get GCD Degree',...
    'Plotting',...
    'Preprocessing',...
    'Results',...
    'Sylvester Matrix');

addpath(genpath('APF'));
addpath(genpath('Examples'));
addpath(genpath('Low Rank Approx'));


% %
% Get Example

[fxy_exact, gxy_exact,dxy_exact,uxy_exact,vxy_exact,m,n,t_exact] = Examples_GCD(ex_num);

[t1,t2] = GetDegree(dxy_exact);
fprintf([mfilename ' : ' sprintf('Total Degree of GCD : %i \n',t_exact)]);
fprintf([mfilename ' : ' sprintf('Deg_x of GCD : %i \n',t1)]);
fprintf([mfilename ' : ' sprintf('Deg_y of GCD : %i \n',t2)]);

% %
% Add Noise to the coefficients

% Add noise to the coefficients of f and g
[fxy_matrix, ~] = AddNoiseToPoly2(fxy_exact,emin,emax);
[gxy_matrix, ~] = AddNoiseToPoly2(gxy_exact,emin,emax);


% %
% %
% Calculate GCD
lower_limit = 1;
upper_limit = min(m,n);

% Calculate the gcd, and quotient polynomials of f(x,y) and g(x,y)
[fxy_calc,gxy_calc,dxy_calc,uxy_calc, vxy_calc,t1,t2] = o_gcd_mymethod(fxy_matrix,gxy_matrix,...
    m,n,[lower_limit,upper_limit]);

% %
% %
% %
% Results.

% Print coefficients of computed d(x,y)
%PrintoutCoefficients('d',dxy_calc,dxy_exact)

% Get error d(x,y)
error.dxy = GetDistance('d',dxy_calc,dxy_exact);
error.uxy = GetDistance('u',uxy_calc,uxy_exact);
error.vxy = GetDistance('v',vxy_calc,vxy_exact);

% Output to file
PrintToFile(m,n,error);


end


function [dist] = GetDistance(name,matrix_calc,matrix_exact)
% Given two matrices, get the distance between them.

% Normalise coefficients 
matrix_calc = normalise(matrix_calc);
matrix_exact = normalise(matrix_exact);

% Get Distance between f(x,y) computed and f(x,y) exact.
dist = (norm(matrix_exact,'fro') - norm(matrix_calc,'fro') )./ norm(matrix_exact,'fro');

% Printout
fprintf([mfilename ' : ' sprintf('Analysis of Coefficients of %s(x,y) computed vs %s(x,y) exact: \n',name,name)])
fprintf([mfilename ' : ' sprintf('Distance between exact and calculated matrix: %2.4e \n',dist)])



end


function []= PrintToFile(m,n,error)

% Global settings
global SETTINGS

% File name
fullFileName = 'Results/Results_o_gcd.txt';


if exist('Results/Results_o_gcd.txt', 'file')
    fileID = fopen('Results/Results_o_gcd.txt','a');
    fprintf(fileID,'%s, \t %s, \t %s,  \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s \n',...
        datetime(),...
        SETTINGS.EX_NUM,...
        num2str(m),...
        num2str(n),...
        error.uxy,...
        error.vxy,...
        error.dxy,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
        SETTINGS.APF_METHOD,...
        SETTINGS.EMIN);
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end

end
