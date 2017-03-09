function [dxy_calc] = o_gcd_Bivariate_2Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_build_method, factorisation_build_method)
% o_gcd_2Polys(ex_num, el, mean_method, bool_alpha_theta, low_rank_approx_method)
%
% Given an example number and set of parameters, obtain GCD of the two
% polynomials f(x,y) and g(x,y) in the given example file.
% Where the polynomials f(x,y) and g(x,y) are defined as polynomails in
% the Bernstein Basis.
%
% % Inputs
%
% ex_num : (String) Example Number
%
% emin : (Float) Lower noise level
%
% emax : (Float) Upper noise level
%
% mean_method (String)
%       'None'
%       'Geometric Mean Matlab Method'
%
% bool_alpha_theta (Boolean)
%       true : Include Preprocessing
%       false : Exclude Preprocessing
%
% low_rank_approx_method (String)
%       'Standard SNTLN' : Include SNTLN
%       'Standard STLN : Include STLN
%       'None' : Exclude SNTLN
%
% apf_method (String)
%       'None'
%       'Standard APF Nonlinear'
%       'Standard APF Linear'
%
% sylvester_build_method (String)
%       'T'
%       'DT'
%       'DTQ'
%       'TQ'
%
% factorisation_build_method (String)
%
%       'HCG'
%       'HC'
%
% % Examples
%
% >> o_gcd_Bivariate_2Polys('1', 1e-12, 1e-10, 'None', false, 'Standard STLN', 'None', 'DTQ', 'HCG')
% >> o_gcd_Bivariate_2Polys('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ','HCG')
% >> o_gcd_Bivariate_2Polys('1', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'None', 'DTQ', 'HCG')

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
    low_rank_approx_method, apf_method, sylvester_build_method, factorisation_build_method)

% Add subfolders
restoredefaultpath

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 

% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

% Print Parameters to console
fprintf('INPUTS. \n')
fprintf('EXAMPLE NUMBER %s \n',ex_num)
fprintf('EMIN : %s \n',emin)
fprintf('EMAX : %s \n',emax)
fprintf('MEAN METHOD : %s \n', mean_method)
fprintf('PREPROCESSING : %s \n',bool_alpha_theta)
fprintf('LOW RANK METHOD : %s \n',low_rank_approx_method)
fprintf('APF METHOD : %s \n', apf_method)
fprintf('Sylvester Matrix Type : %s \n', sylvester_build_method);
fprintf('Factorisation Matrix Type : %s \n', factorisation_build_method);

% Get Example polynomials
[fxy_exact, gxy_exact, dxy_exact, uxy_exact, vxy_exact, m, n, t_exact] = Examples_GCD_Bivariate_2Polys(ex_num);

% Get degree of f(x,y), g(x,y) and d(x,y)
[t1, t2] = GetDegree_Bivariate(dxy_exact);
[m1, m2] = GetDegree_Bivariate(fxy_exact);
[n1, n2] = GetDegree_Bivariate(gxy_exact);

% Print degree structure
fprintf('The Degree Structure of f(x) : m = %i \t m1 = %i \t m2 = %i \n', m, m1, m2)
fprintf('The Degree Structure of g(x) : n = %i \t n1 = %i \t n2 = %i \n', n, n1, n2)
fprintf('The Degree Structure of d(x) : t = %i \t t1 = %i \t t2 = %i \n', t_exact, t1, t2)

% Add noise to the coefficients of f(x,y) and g(x,y)
[fxy_matrix, ~] = AddVariableNoiseToPoly(fxy_exact, emin, emax);
[gxy_matrix, ~] = AddVariableNoiseToPoly(gxy_exact, emin, emax);

% Calculate GCD

% Set upper and lower limit for the degree structure of d(x,y)
lowerLimit_t1 = 0;
upperLimit_t1 = min(m1, n1);
lowerLimit_t2 = 0;
upperLimit_t2 = min(m2, n2);
limits_t1 = [lowerLimit_t1 upperLimit_t1];
limits_t2 = [lowerLimit_t2 upperLimit_t2];

% Calculate the gcd, and quotient polynomials of f(x,y) and g(x,y)
[~, ~, dxy_calc, uxy_calc, vxy_calc, t1, t2] = ...
    o_gcd_mymethod_Bivariate_2Polys(fxy_matrix, gxy_matrix, limits_t1, limits_t2);

% Get error in coefficients of d(x,y), u(x,y) and v(x,y)
myError.dxy = GetDistance('d', dxy_calc, dxy_exact);
myError.uxy = GetDistance('u', uxy_calc, uxy_exact);
myError.vxy = GetDistance('v', vxy_calc, vxy_exact);
 
% Output results to file
PrintToFile(m, n, t1, t2, myError);


end


function [dist] = GetDistance(name, fxy_calc, fxy_exact)
% Given two matrices, get the distance between them.
%
% name : (String)
%
% fxy_calc : (Matrix) Coefficients of f(x,y) as calculated
%
% fxy_exact : (Matrix) Coefficients of f(x,y) in exact form

% Normalise coefficients
fxy_calc = normalise(fxy_calc);
fxy_exact = normalise(fxy_exact);

try
    % Get Distance between f(x,y) computed and f(x,y) exact.
    dist = norm(fxy_exact-fxy_calc,'fro') ./ norm(fxy_exact,'fro');

catch
    
    dist = 1000;
end

% Printout
fprintf([mfilename ' : ' sprintf('Analysis of Coefficients of %s(x,y) computed vs %s(x,y) exact: \n',name,name)])
fprintf([mfilename ' : ' sprintf('Distance between exact and calculated matrix: %2.4e \n',dist)])



end


function []= PrintToFile(m, n, t1, t2, myError)
% Print the results to a file
%
% % Inputs
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomail g(x,y)
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
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime(),...
            SETTINGS.EX_NUM,...
            num2str(m),...
            num2str(n),...
            num2str(t1),...
            num2str(t2),...
            num2str(myError.uxy),...
            num2str(myError.vxy),...
            num2str(myError.dxy),...
            SETTINGS.MEAN_METHOD,...
            num2str(SETTINGS.BOOL_ALPHA_THETA),...
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
        fprintf(fileID,'DATE,EX_NUM,m,n,t1,t2,ERROR_UX,ERROR_VX,ERROR_DX,MEAN_METHOD,BOOL_ALPHA_THETA, EMIN, EMAX, LOW_RANK_APPROX_METHOD,LOW_RANK_ITE, APF_METHOD, APF_ITE,sylvester_build_method \n');
    end



end
