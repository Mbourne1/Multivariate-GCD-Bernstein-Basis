function [dxy_calc] = o_gcd_Bivariate_3Polys(ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, factorisation_build_method, ...
    rank_revealing_metric, nEquations, degree_method)
% o_gcd(ex_num, el, mean_method, bool_alpha_theta, ...
% low_rank_approx_method, apf_method, sylvester_build_method, ...
% factorisation_build_method, rank_revealing_metric)
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
% mean_method : (String)
%   'None'
%   'Geometric Mean Matlab Method'
%
% bool_alpha_theta : (Boolean)
%   true : Include Preprocessing
%   false : Exclude Preprocessing
%
% low_rank_approx_method : (String)
%   'Standard SNTLN' : Include SNTLN
%   'Standard STLN : Include STLN
%   'None' : Exclude SNTLN
%
% apf_method : (String)
%   'None'
%   'Standard APF Nonlinear'
%   'Standard APF Linear'
%
% sylvester_matrix_variant : (String)
%   'T'
%   'DT'
%   'DTQ'
%   'TQ'
%
% factorisation_build_method : (String)
%   'HCG' : 
%   'HC' :
%
% rank_revealing_metric : (String)
%   'Minimum Singular Values'
%
%
%
%
% % Examples
%
% >> o_gcd_Bivariate_3Polys('1',1e-12,1e-10,'None',false,'None','None','DTQ', 'HCG', 'Minimum Singular Values')
% >> o_gcd_Bivariate_3Polys('1',1e-12,1e-10,'Geometric Mean Matlab Method',true,'Standard STLN','None','DTQ', 'HCG', 'Minimum Singular Values')




% %
% Set Variables

% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end


SetGlobalVariables_GCD_3Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_matrix_variant, factorisation_build_method, rank_revealing_metric, nEquations);


% Print Parameters to console
fprintf('INPUTS. \n')
fprintf('EXAMPLE NUMBER %s \n',ex_num)
fprintf('EMIN : %s \n',emin)
fprintf('EMAX : %s \n',emax)
fprintf('MEAN METHOD : %s \n', mean_method)
fprintf('PREPROCESSING : %s \n', num2str(bool_alpha_theta))
fprintf('LOW RANK METHOD : %s \n',low_rank_approx_method)
fprintf('APF METHOD : %s \n', apf_method)
fprintf('SYLVESTER MATRIX VARIANT : %s \n', sylvester_matrix_variant)
fprintf('RANK REVEALING METRIC : %s \n', rank_revealing_metric)

% %
% Get Example
[fxy_exact, gxy_exact, hxy_exact, dxy_exact, uxy_exact, vxy_exact, ...
    wxy_exact, m, n, o, t_exact] = Examples_GCD_Bivariate_3Polys(ex_num);

% Get the degree of the GCD
[t1, t2] = GetDegree_Bivariate(dxy_exact);

fprintf([mfilename ' : ' sprintf('Total Degree of GCD : %i \n',t_exact)]);
fprintf([mfilename ' : ' sprintf('Deg_x of GCD : %i \n',t1)]);
fprintf([mfilename ' : ' sprintf('Deg_y of GCD : %i \n',t2)]);

% %
% Add Noise to the coefficients

% Add noise to the coefficients of f and g
[fxy, ~] = AddVariableNoiseToPoly(fxy_exact, emin, emax);
[gxy, ~] = AddVariableNoiseToPoly(gxy_exact, emin, emax);
[hxy, ~] = AddVariableNoiseToPoly(hxy_exact, emin, emax);


% Get degree of f(x,y), g(x,y) and h(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);
[t1_exact, t2_exact] = GetDegree_Bivariate(dxy_exact);


% Calculate GCD
lowerLimit_t1 = 1;
upperLimit_t1 = min([m1, n1, o1]);
limits_t1 = [lowerLimit_t1, upperLimit_t1];

lowerLimit_t2 = 1;
upperLimit_t2 = min([m2, n2, o2]);
limits_t2 = [lowerLimit_t2, upperLimit_t2];


rank_range = [0,0];

% Calculate the gcd, and quotient polynomials of f(x,y) and g(x,y)
[fxy_calc, gxy_calc, hxy_exact, dxy_calc, uxy_calc, vxy_calc, wxy_calc, ...
    t1, t2] = o_gcd_mymethod_Bivariate_3Polys(fxy, gxy, hxy, limits_t1, ...
    limits_t2, rank_range, degree_method);

% %
% %
% %
% Results.

% Print coefficients of computed d(x,y)
%PrintoutCoefficients('d',dxy_calc,dxy_exact)

% Get error d(x,y)
error.dxy = GetDistance('d', dxy_calc, dxy_exact);
error.uxy = GetDistance('u', uxy_calc, uxy_exact);
error.vxy = GetDistance('v', vxy_calc, vxy_exact);
error.wxy = GetDistance('w', wxy_calc, wxy_exact);

% Output to file
PrintToFile(m, n, o, t1, t2, t1_exact, t2_exact, error);


end


function [dist] = GetDistance(name, matrix_calc, matrix_exact)
% Given two matrices, get the distance between them.
%
% % Inputs
%
% name : (String)
% 
% matrix_calc : (Matrix) The inexact matrix
%
% matrix_exact : (Matrix) The exact matrix
%
% % Outputs
%
% dist : (Float) Distance between the exact and inexact matrix


% Normalise coefficients
matrix_calc = normalise(matrix_calc);
matrix_exact = normalise(matrix_exact);

try
    % Get Distance between f(x,y) computed and f(x,y) exact.
    dist = norm(matrix_exact - matrix_calc,'fro') ./ norm(matrix_exact,'fro');
catch
    dist = 1000;
end

% Printout
fprintf([mfilename ' : ' sprintf('Analysis of Coefficients of %s(x,y) computed vs %s(x,y) exact: \n',name,name)])
fprintf([mfilename ' : ' sprintf('Distance between exact and calculated matrix: %2.4e \n',dist)])



end


function []= PrintToFile(m, n, o, t1, t2, t1_exact, t2_exact, error)
% Print the results to a file
%
% % Inputs
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomail g(x,y)
%
% t1 : (Int) The computed Degree of d(x,y) with respect to x
%
% t2  : (Int) The comuted degree of d(x,y) with respect to y
%
% t1_exact : (Int)
% error: (Float Float Float Float) Contains error.uxy, error.vxy error.wxy and error.dxy




% Global settings
global SETTINGS

fullFileName = sprintf('Results/Results_o_gcd_3Polys.dat');

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
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime(),...
            SETTINGS.EX_NUM,...
            num2str(m),...
            num2str(n),...
            num2str(o),...
            num2str(t1),...
            num2str(t2),...
            num2str(t1_exact),...
            num2str(t2_exact),...
            num2str(error.uxy),...
            num2str(error.vxy),...
            num2str(error.wxy),...
            num2str(error.dxy),...
            SETTINGS.MEAN_METHOD,...
            num2str(SETTINGS.BOOL_ALPHA_THETA),...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
            num2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
            SETTINGS.APF_METHOD,...
            num2str(SETTINGS.APF_REQ_ITE),...
            SETTINGS.SYLVESTER_MATRIX_VARIANT...
            );
    end

    function WriteHeader()
        
        strHeader = ['DATE, EX_NUM, m, n, o, t1, t2, t1_exact, t2_exact,'...
            'ERROR_UXY, ERROR_VXY, ERROR_WXY, ERROR_DX, '...
            'MEAN_METHOD, BOOL_ALPHA_THETA, EMIN, EMAX, '...
            'LOW_RANK_APPROX_METHOD, LOW_RANK_ITE, APF_METHOD, '...
            'APF_ITE, SYLVESTER_MATRIX_VARIANT \n'];
        
        fprintf(fileID, strHeader);
    end



end
