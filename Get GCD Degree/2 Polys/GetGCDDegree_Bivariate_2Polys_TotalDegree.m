function [t1, t2, GM_fx, GM_gx, alpha, th1, th2] = GetGCDDegree_Bivariate_2Polys_TotalDegree(fxy, gxy, myLimits_t1, myLimits_t2, limits_t1, limits_t2)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y). This is done by first setting k_{1} = 1.
% Computing the set of Sylvester subresultant matrices S_{k1,k2} for all
% possible k2 values, from which t_{2} (the degree of the gcd with respect
% to y) is computed. Then construct the set of Sylvester subresultant
% matrices S_{k1,t_{2}} for all possible k_{1} values, to determine the
% degree of the GCD with respect to x
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) in the Bernstein basis
%
% gxy : (Matrix) Coefficients of polynomial g(x,y) in the Bernstein basis
%
% myLimits_t1 : [(Int) (Int)] 
%
% myLimits_t2 : [(Int) (Int)]
%
% limits_t1 : [(Int) (Int)] [lowerLimit upperLimit] when computing the
% degree of the GCD with respect to x
%
% limits_t2 : [(Int) (int)] [lowerLimit upperLimit] when computing the
% degree of the GCD with respect to y
%
% % Outputs
%
% t1 : (Int) Degree of d(x,y) with respect to x
%
% t2 : (Int) Degree of d(x,y) with respect to y
%
% GM_fx : (Float) Geometric mean of entries in first partition of Syvlester
%   matrix
%
% GM_gx : (Float) Geometric mean of Entries in second partition of Sylvester matrix
%
% alpha : (Float) Optimal value of alpha
%
% th1 : (Float) Optimal value of theta_{1}
%
% th2 : (Float) Optimal value of theta_{2}



% Get my lower limit for the computation of t_{1} and t_{2}
myLowerLimit_t1 = myLimits_t1(1);
myLowerLimit_t2 = myLimits_t2(1);


% Set x = 1, so
x = 1;

% Build the set of Sylvester subresultant matrices S_{1,k2} for all k_{2}
% values.
[arr_R1, arr_SingularValues] = getSubresultants(fxy, gxy, myLimits_t2, 'x', x);

% Return a metric used to determine the degree of the GCD with respect to y
[metric] = getMetric(arr_R1, arr_SingularValues, myLimits_t2);

% Get degree of GCD with respect to y
vDelta_y = diff(log10(metric));
[~, index] = max(vDelta_y);
t2 = index + myLowerLimit_t2 -1;


% Build the set of Sylvester subresultant matrices S_{k1, t_{2}}(f,g)
[arr_R1, arr_SingularValues, vAlpha, vTheta1, vTheta2, vGM_fx, vGM_gx] = getSubresultants(fxy, gxy, myLimits_t1, 'y', t2);

% Return a metric used to determin the degree of the GCD with respect to x
[metric] = getMetric(arr_R1, arr_SingularValues, myLimits_t1, myLimits_t2, limits_t1, limits_t2);

% Get the degree of the GCD with respect to x
vDelta_x = diff(log10(metric));
[~, index] = max(vDelta_x);
t1 = index + myLowerLimit_t1 -1;


% Outputs
GM_fx = vGM_fx(t1 - myLowerLimit_t1 +1);
GM_gx = vGM_gx(t1 - myLowerLimit_t1 +1);
alpha = vAlpha(t1 - myLowerLimit_t1 +1);
th1 = vTheta1(t1 - myLowerLimit_t1 +1);
th2 = vTheta2(t1 - myLowerLimit_t1 +1);



% Print results to console.
LineBreakMedium()
fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n',t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n',t2)])
LineBreakMedium()

end


function [arr_R1, arr_SingularValues, vAlpha, vTheta1, vTheta2, vGM_fx, vGM_gx] = getSubresultants(fxy, gxy, limits_k, str_fixed_var, fixed_var)
% Compute the degree (Either t_{1} or t_{2}) of the GCD d(x,y) given that
% by constructing the sequence of Sylvester subresultant matrices S_{k1,k2}
% where either k_{1} or k_{2} is constant.
%
% % Inputs
%
% fxy : (Matrix) : Coefficients of f(x,y)
%
% gxy : (Matrix) : Coefficients of g(x,y)
%
% limits_k : [Int) (Int)]
%
% str_fixed_var : (String) Either 'x' or 'y'
%
% str_var :  (Integer) Value of fixed variable.

% Get lower and upper limit in computation of the degree of the GCD.
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% Get number of Sylvester subresultant matrices to be constructed
nSubresultants = upperLimit_k - lowerLimit_k + 1;


% Initialise an array to store singular values of each of the Sylvester
% subresultant matrices, and an array to store matrix R1 from QR
% decompostiion of S_{k1,k2}
% matrix R1
arr_SingularValues = cell(nSubresultants, 1);
arr_R1 = cell(nSubresultants, 1);

% Initialise some vectors
vAlpha = zeros(nSubresultants, 1);
vTheta1 = zeros(nSubresultants, 1);
vTheta2 = zeros(nSubresultants, 1);
vGM_fx = zeros(nSubresultants, 1);
vGM_gx = zeros(nSubresultants, 1);



% For each Sylvester subresultant matrix
for i2 = 1:1:nSubresultants
    
    switch str_fixed_var
        case 'x'
            
            k1 = fixed_var;
            k2 = lowerLimit_k + (i2-1);
            
        case 'y'
            
            k1 = (i2-1);
            k2 = fixed_var;
            
    end
    
    
    % Preprocess the matrix S_{k1,k2}
    [GM_fx, GM_gx, alpha, th1, th2] = Preprocess(fxy, gxy, k1, k2);
    
    % Divide f(x) by geometric mean
    fxy_matrix_n = fxy ./ GM_fx;
    
    % Divide g(x) by geometric mean
    gxy_matrix_n = gxy ./ GM_gx;
    
    % Get f(w,w) from f(x,y)
    fww_matrix = GetWithThetas(fxy_matrix_n, th1, th2);
    
    % Get g(w,w) from g(x,y)
    gww_matrix = GetWithThetas(gxy_matrix_n, th1, th2);
    
    % Build the Sylvester subresultant matrix S_{k1,k2}(f,g)
    Sk1k2 = BuildSubresultant_Bivariate_2Polys(fww_matrix, alpha.*gww_matrix, k1, k2);
    
    % Get the singular values of S_{k1,k2}
    arr_SingularValues{i2} = svd(Sk1k2);
    
    % Get QR decomposition of S_{k1,k2}
    [~,R] = qr(Sk1k2);
    [~,c] = size(R);
    arr_R1{i2} = R(1:c,1:c);
    
    % Store the optimal alpha
    vAlpha(i2) = alpha;
    
    % Store the optimal theta_{1}
    vTheta1(i2) = th1;
    
    % Store the optimal theta_{2}
    vTheta2(i2) = th2;
    
    % Store the optimal lambda
    vGM_fx(i2) = GM_fx;
    
    % Store the optimal mu
    vGM_gx(i2) = GM_gx;
    
end



end

function [metric] = getMetric(arr_R1, arr_SingularValues, limits_k1, limits_t1)
% This function returns a metric from which the degree of the GCD can be
% computed (Either with respect to x or y). Note the remainder of this
% function is worded to suggest we compute the degree of the GCD with
% respect to x, but can equally be applied to computing the degree with
% respect to y.
%
% % Inputs
%
% arr_R1 : (1D Array of Matrices) : Contains array of matrices R_{1} from
% QR decomposition of Sylvester subresultant matrices S_{k1,k2}, where
% either k1 or k2 is constant.
%
% arr_SingularValues : (1D Array of Vectors) : Vectors containing Singular
% values of each Sylvester subresultant matrix S_{k1,k2} where either k_{1}
% or k_{2} is constant
%
% limits_k1 : (Int Int) Set limits for the computation of the degree of the GCD,
%
%
% limits_t1 :
%
% % Outputs
%
% metric (Vector) Contains a metric from which the degree of the GCD is
% determined by the max change in the vector.

% Get upper and lower limit of
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);

% Get number of Sylvester subresultant matrices required to compute degree
% of GCD.
nSubresultants = upperLimit_k1 - lowerLimit_k1 + 1;

global SETTINGS
switch SETTINGS.RANK_REVEALING_METRIC
    
    case 'R1 Row Norms'
        
        % Initialise matrices to store max and minimum row norms
        arr_R1_RowNorms = cell(nSubresultants, 1);
        vMaxRowNorm = zeros(nSubresultants, 1);
        vMinRowNorm = zeros(nSubresultants, 1);
        
        % Get maximum and minimum row norms of R_{1} from QR decompositon
        % of S_{k1,k2}
        
        for i1 = 1:1:nSubresultants
            
            
            arr_R1_RowNorms{i1} = sqrt(sum(arr_R1{i1}.^2,2))./norm(arr_R1{i1});
            vMaxRowNorm(i1) = max(arr_R1_RowNorms{i1});
            vMinRowNorm(i1) = min(arr_R1_RowNorms{i1});
            
        end
        
        % Plot graphs
        if(SETTINGS.PLOT_GRAPHS)
            
            % Plot graphs
            %plotR1RowNorms_1Dimensional(arr_R1_RowNorms, myLimits_t1, limits_t1)
            plotMaxMinRowNorms_1Dimensional(vMaxRowNorm, vMinRowNorm, limits_k1, limits_t1);
            
        end
        
        % Set rank revealing metric
        metric = vMinRowNorm ./ vMaxRowNorm;
         
    case 'R1 Row Diagonals'
        
        % Initialise matrices to store max and minimum diagonals
        vMaxDiagonal_R1 = zeros(nSubresultants, 1);
        vMinDiagonal_R1 = zeros(nSubresultants, 1);
        
        % Get maximum and minimum diagonal of each R_{k1,k2} from QR
        % decomposition of S_{k1,k2}
        
        
        for i1 = 1:1:nSubresultants
            
            vMaxDiagonal_R1(i1) = max(abs(diag(arr_R1{i1})));
            vMinDiagonal_R1(i1) = min(abs(diag(arr_R1{i1})));
            
        end
        
        
        % Plot graphs
        if(SETTINGS.PLOT_GRAPHS)
            
            %plotR1Diagonals_1Dimensional(arr_R1, myLimits_t1, limits_t1);
            plotMaxMinDiagonalsR1_1Dimensional(vMaxDiagonal_R1, vMinDiagonal_R1, limits_k1, limits_t1);
        end
        
        % Set rank revealing metric
        metric = vMinDiagonal_R1 ./ vMaxDiagonal_R1;
        
    case 'Singular Values'
        
        % Initialise matrix to store minimum singular values of SVD of each
        % S_{k1,k2}
        vMinimumSingularValues = zeros(nSubresultants,1);
        
        for i1 = 1:1:nSubresultants
            
            vMinimumSingularValues(i1) = min(arr_SingularValues{i1});
            
        end
        
        
        % Plot graphs
        if(SETTINGS.PLOT_GRAPHS)
            
            %plotSingularValues_1Dimensional(arr_SingularValues, myLimits_t1, limits_t1)
            plotMinimumSingularValues_1Dimensional(vMinimumSingularValues, limits_k1, limits_t1)
            
        end
        
        
        % Set rank revealing metric
        metric = vMinimumSingularValues;
        
    case 'Residuals'
        
        error('Error : Code not yet developed for this branch')
        
    otherwise
        
        error('Error : Not a valid metric');
        
end


end

