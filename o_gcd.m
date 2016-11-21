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
% emin - Lower noise level
%
% emax - Upper noise level
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

% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

% Set global variables
SetGlobalVariables(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method)

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

% Print Parameters to console
fprintf('INPUTS. \n')
fprintf('EXAMPLE NUMBER %s \n',ex_num)
fprintf('EMIN : %s \n',emin)
fprintf('EMAX : %s \n',emax)
fprintf('MEAN METHOD : %s \n', mean_method)
fprintf('PREPROCESSING : %s \n',bool_alpha_theta)
fprintf('LOW RANK METHOD : %s \n',low_rank_approx_method)
fprintf('APF METHOD : %s \n', apf_method)

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
[fxy_matrix, ~] = AddNoiseToPoly(fxy_exact,emin);
[gxy_matrix, ~] = AddNoiseToPoly(gxy_exact,emin);


% %
% %
% Calculate GCD
lower_limit = 1;
upper_limit = min(m,n);

% Calculate the gcd, and quotient polynomials of f(x,y) and g(x,y)
[fxy_calc,gxy_calc,dxy_calc,uxy_calc, vxy_calc, t1, t2] = ...
    o_gcd_mymethod(fxy_matrix, gxy_matrix, m, n, [lower_limit,upper_limit]);

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

try
    % Get Distance between f(x,y) computed and f(x,y) exact.
    dist = (norm(matrix_exact,'fro') - norm(matrix_calc,'fro') )./ norm(matrix_exact,'fro');
catch
    dist = 1000;
end

% Printout
fprintf([mfilename ' : ' sprintf('Analysis of Coefficients of %s(x,y) computed vs %s(x,y) exact: \n',name,name)])
fprintf([mfilename ' : ' sprintf('Distance between exact and calculated matrix: %2.4e \n',dist)])



end


function []= PrintToFile(m,n,error)
% Print the results to a file
%
% % Inputs
%
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomail g(x,y)
%
% error: Contains error.uxy, error.vxy and error.dxy


% Global settings
global SETTINGS

fullFileName = sprintf('Results/Results_o_gcd%s.txt',datetime('today'));

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    WriteNewLine()
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    WriteNewLine()
    fclose(fileID);
    
end

    function WriteNewLine()
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime(),...
            SETTINGS.EX_NUM,...
            num2str(m),...
            num2str(n),...
            num2str(error.uxy),...
            num2str(error.vxy),...
            num2str(error.dxy),...
            SETTINGS.MEAN_METHOD,...
            SETTINGS.BOOL_ALPHA_THETA,...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
            num2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
            SETTINGS.APF_METHOD,...
            num2str(SETTINGS.APF_REQ_ITE),...
            SETTINGS.SYLVESTER_BUILD_METHOD...
            );
    end

    function WriteHeader()
        fprintf(fileID,'DATE,EX_NUM,m,n,ERROR_UX,ERROR_VX,ERROR_DX,MEAN_METHOD,BOOL_ALPHA_THETA, EMIN, EMAX, LOW_RANK_APPROX_METHOD,LOW_RANK_ITE, APF_METHOD, APF_ITE,sylvester_build_method \n');
    end



end
