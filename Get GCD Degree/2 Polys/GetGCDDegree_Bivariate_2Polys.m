function [t1, t2, GM_fx, GM_gx, alpha, th1, th2] = GetGCDDegree_Bivariate_2Polys(fxy, gxy, limits_t1, limits_t2)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y)
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) in the Bernstein basis
%
% gxy : (Matrix) Coefficients of polynomial g(x,y) in the Bernstein basis
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


%
global SETTINGS


% Get the degree structure of polynomial f(x,y) and g(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Set my limits for the computation of the degree of the GCD.
limits_k1 = [0 min(m1,n1)];
limits_k2 = [0 min(m2,n2)];

lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);
lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

nSubresultants_k1 = upperLimit_k1 - lowerLimit_k1 + 1;
nSubresultants_k2 = upperLimit_k2 - lowerLimit_k2 + 1;

% Initialise some vectors
arr_SingularValues = cell(nSubresultants_k1, nSubresultants_k2);


matrix_Alpha = zeros(nSubresultants_k1, nSubresultants_k2);
matrix_Theta1 = zeros(nSubresultants_k1, nSubresultants_k2);
matrix_Theta2 = zeros(nSubresultants_k1, nSubresultants_k2);
matrix_GM_fx = zeros(nSubresultants_k1, nSubresultants_k2);
matrix_GM_gx = zeros(nSubresultants_k1, nSubresultants_k2);

arr_R1 = cell(nSubresultants_k1, nSubresultants_k2);


% For every row in the matrix
for i1 = 1:1:nSubresultants_k1
    
    for i2 = 1:1:nSubresultants_k2
        
        k1 = lowerLimit_k1 + (i1-1);
        k2 = lowerLimit_k2 + (i2-1);
        
        % Preprocessing
        [GM_fx, GM_gx, alpha, th1, th2] = Preprocess_Bivariate_2Polys(fxy, gxy, k1, k2);
        
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
        
        % Get the singular values
        arr_SingularValues{i1, i2} = svd(Sk1k2);
        
        % Get QR decomposition
        [~,R] = qr(Sk1k2);
        [~,c] = size(R);
        arr_R1{i1, i2} = R(1:c,1:c);
        
        
        % Store the optimal alpha
        matrix_Alpha(i1, i2) = alpha;
        
        % Store the optimal theta_{1}
        matrix_Theta1(i1, i2) = th1;
        
        % Store the optimal theta_{2}
        matrix_Theta2(i1, i2) = th2;
        
        % Store the optimal lambda
        matrix_GM_fx(i1, i2) = GM_fx;
        
        % Store the optimal mu
        matrix_GM_gx(i1, i2) = GM_gx;
        
    end
end

% if (SETTINGS.PLOT_GRAPHS)
%     
%     x_vec = lowerLimit_k1:1:upperLimit_k1;
%     y_vec = lowerLimit_k2:1:upperLimit_k2;
%     
%     [X,Y] = meshgrid(x_vec,y_vec);
%     
%     figure_name = sprintf('Geometric Mean of f(x) %s', SETTINGS.SYLVESTER_BUILD_METHOD);
%     figure('name',figure_name)
%     hold on
%     mesh(X,Y,log10(matrix_GM_fx)');
%     hold off
%     
%     figure_name = sprintf('Geometric Mean of g(x) in %s', SETTINGS.SYLVESTER_BUILD_METHOD);
%     figure('name', figure_name)
%     hold on
%     mesh(X,Y,log10(matrix_GM_gx)');
%     hold off
%     
%     figure_name = sprintf('Optimal values of alpha, th2 and th2 in %s', SETTINGS.SYLVESTER_BUILD_METHOD);
%     figure('name',figure_name)
%     hold on
%     %mesh(X, Y, log10(matrix_Alpha)');
%     mesh(X, Y, log10(matrix_Theta1)');
%     mesh(X, Y, log10(matrix_Theta2)');
%     hold off
% end



% Metric used to compute the degree of the GCD
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals

switch SETTINGS.RANK_REVEALING_METRIC
    case 'R1 Row Norms'
        
        % Initialise matrices to store max and minimum row norms
        arr_R1_RowNorms = cell(nSubresultants_k1, nSubresultants_k2);
        mat_MaxRowNorm = zeros(nSubresultants_k1, nSubresultants_k2);
        mat_MinRowNorm = zeros(nSubresultants_k1, nSubresultants_k2);
        
        % Get maximum and minimum row norms of R_{1} from QR decompositon
        % of S_{k1,k2}
        for i1 = 1:1:nSubresultants_k1
            
            for i2 = 1:1:nSubresultants_k2
                
                
                arr_R1_RowNorms{i1, i2} = sqrt(sum(arr_R1{i1, i2}.^2,2))./norm(arr_R1{i1, i2});
                mat_MaxRowNorm(i1, i2) = max(arr_R1_RowNorms{i1, i2});
                mat_MinRowNorm(i1, i2) = min(arr_R1_RowNorms{i1, i2});
                
            end
        end
        
        
        % Plot graphs
        if(SETTINGS.PLOT_GRAPHS)
        
            plotR1RowNorms(arr_R1_RowNorms, limits_k1, limits_k2, limits_t1, limits_t2)
            plotMaxMinRowNorms(mat_MaxRowNorm, mat_MinRowNorm, limits_k1, limits_k2, limits_t1, limits_t2);
        
        end
        
        % Set Metric
        metric = mat_MinRowNorm./mat_MaxRowNorm;
        
    case 'R1 Row Diagonals'
        
        % Initialise matrices to store max and minimum diagonals
        mat_MaxDiagonal_R1 = zeros(nSubresultants_k1, nSubresultants_k2);
        mat_MinDiagonal_R1 = zeros(nSubresultants_k1, nSubresultants_k2);
        
        % Get maximum and minimum diagonal of each R_{k1,k2} from QR
        % decomposition of S_{k1,k2}
        for i1 = 1:1:nSubresultants_k1
            
            for i2 = 1:1:nSubresultants_k2
                
                mat_MaxDiagonal_R1(i1,i2) = max(abs(diag(arr_R1{i1,i2})));
                mat_MinDiagonal_R1(i1,i2) = min(abs(diag(arr_R1{i1,i2})));
                
            end
        end
        
        % Plot graphs
        if(SETTINGS.PLOT_GRAPHS)
        
            plotR1Diagonals(arr_R1, limits_k1, limits_k2, limits_t1, limits_t2);
            plotMaxMinDiagonalsR1(mat_MaxDiagonal_R1, mat_MinDiagonal_R1, limits_k1, limits_k2, limits_t1, limits_t2);
        
        end
        
        % Set metric
        metric = mat_MinDiagonal_R1 ./ mat_MaxDiagonal_R1;
        
    case 'Singular Values'
        
        % Initialise matrix to store minimum singular values of SVD of each
        % S_{k1,k2}
        mat_MinimumSingularValues = zeros(min(m1,n1)+1, min(m2,n2)+1);
        
        for i1 = 1:1:nSubresultants_k1
            
            for i2 = 1:1:nSubresultants_k2
                
                
                % k1 = lowerLimit_k1 + (i1 - 1)
                % k2 = lowerLimit_k2 + (i2 - 1)
                mat_MinimumSingularValues(i1,i2) = min(min(arr_SingularValues{i1,i2}));
                
            end
        end
        
        % Plot Graphs
        if(SETTINGS.PLOT_GRAPHS)
            plotSingularValues(arr_SingularValues, limits_k1, limits_k2, limits_t1, limits_t2)
            plotMinimumSingularValues(mat_MinimumSingularValues, limits_k1, limits_k2, limits_t1, limits_t2)
        end
        
        
        
        % Set metric
        metric = mat_MinimumSingularValues;
        
    case 'Residuals'
        
        error('Error : Code not yet developed for this branch')
        
    otherwise
        
        error('Error : Not a valid metric');
        
end



% Compute the degree of the GCD
delta_x = diff(log10(metric),1,1);
vDelta_x = sum(delta_x,2);
[~, index] = max(vDelta_x);
t1 = index + lowerLimit_k1 -1;

delta_y = diff(log10(metric),1,2);
vDelta_y = sum(delta_y,1);
[~, index] = max(vDelta_y);
t2 = index + lowerLimit_k1 -1;


% Outputs
GM_fx = matrix_GM_fx(t1 - lowerLimit_k1 +1, t2 - lowerLimit_k2 +1);
GM_gx = matrix_GM_gx(t1 - lowerLimit_k1 +1, t2 - lowerLimit_k2 +1);
alpha = matrix_Alpha(t1 - lowerLimit_k1 +1, t2 - lowerLimit_k2 +1);
th1 = matrix_Theta1(t1 - lowerLimit_k1 +1, t2 - lowerLimit_k2 +1);
th2 = matrix_Theta2(t1 - lowerLimit_k1 +1, t2 - lowerLimit_k2 +1);

LineBreakMedium()
fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n',t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n',t2)])
LineBreakMedium()

end




