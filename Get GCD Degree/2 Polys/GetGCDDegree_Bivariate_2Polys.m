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


% Get number of subresultant matrices
nSubresultants_k1 = upperLimit_k1 - lowerLimit_k1 + 1;
nSubresultants_k2 = upperLimit_k2 - lowerLimit_k2 + 1;

% Initialise matrices to store values
matrix_Alpha = zeros(nSubresultants_k1, nSubresultants_k2);
matrix_Theta1 = zeros(nSubresultants_k1, nSubresultants_k2);
matrix_Theta2 = zeros(nSubresultants_k1, nSubresultants_k2);
matrix_GM_fx = zeros(nSubresultants_k1, nSubresultants_k2);
matrix_GM_gx = zeros(nSubresultants_k1, nSubresultants_k2);

arr_Sk1k2 = cell(nSubresultants_k1, nSubresultants_k2);
arr_R1 = cell(nSubresultants_k1, nSubresultants_k2);


% For every row in the matrix
for i1 = 1 : 1 : nSubresultants_k1
    
    for i2 = 1 : 1 : nSubresultants_k2
        
        k1 = lowerLimit_k1 + (i1 - 1);
        k2 = lowerLimit_k2 + (i2 - 1);
        
        % Get geometric means of entries of f(x,y) and g(x,y) in the
        % k_{1},k_{2} -th subresultant matrix, and get optimal values of
        % alpha, theta_{1} and theta_{2}
        [GM_fx, GM_gx, alpha, th1, th2] = Preprocess_Bivariate_2Polys(fxy, gxy, k1, k2);
        
        % Divide f(x,y) and g(x,y) by geometric mean
        fxy_matrix_n = fxy ./ GM_fx;
        gxy_matrix_n = gxy ./ GM_gx;
        
        % Get f(w,w) from f(x,y) and g(w,w) from g(x,y)
        fww_matrix = GetWithThetas(fxy_matrix_n, th1, th2);
        gww_matrix = GetWithThetas(gxy_matrix_n, th1, th2);
        
        % Build the subresultant matrix S_{k1, k2}(f,g)
        arr_Sk1k2{i1, i2} = BuildSubresultant_Bivariate_2Polys(fww_matrix, ...
            alpha .* gww_matrix, k1, k2);
        
        % Store the optimal alpha, th1, th2, gm_fx, gm_gx
        matrix_Alpha(i1, i2) = alpha;
        matrix_Theta1(i1, i2) = th1;
        matrix_Theta2(i1, i2) = th2;
        matrix_GM_fx(i1, i2) = GM_fx;
        matrix_GM_gx(i1, i2) = GM_gx;
        
    end
end


if SETTINGS.PLOT_GRAPHS_PREPROCESSING == true
    PlotAlphas(matrix_Alpha, limits_k1, limits_k2, limits_t1, limits_t2);
end


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
                
                [~,R] = qr(arr_Sk1k2{i1, i2});
                [~,c] = size(R);
                arr_R1{i1, i2} = R(1 : c, 1 : c);
                
                
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
        metric = mat_MinRowNorm ./ mat_MaxRowNorm;
        
    case 'R1 Row Diagonals'
        
        % Initialise matrices to store max and minimum diagonals
        mat_MaxDiagonal_R1 = zeros(nSubresultants_k1, nSubresultants_k2);
        mat_MinDiagonal_R1 = zeros(nSubresultants_k1, nSubresultants_k2);
        
        % Get maximum and minimum diagonal of each R_{k1,k2} from QR
        % decomposition of S_{k1,k2}
        for i1 = 1:1:nSubresultants_k1
            
            for i2 = 1:1:nSubresultants_k2
                
                [~,R] = qr(arr_Sk1k2{i1,i2});
                [~,c] = size(R);
                arr_R1{i1, i2} = R(1:c,1:c);
                
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
        
    case 'Minimum Singular Values'
        
        
        
        % Initialise some vectors
        arr_SingularValues = cell(nSubresultants_k1, nSubresultants_k2);
        
        % Initialise matrix to store minimum singular values of SVD of each
        % S_{k1,k2}
        mat_MinimumSingularValues = zeros(min(m1,n1)+1, min(m2,n2)+1);
        
        
        for i1 = 1 : 1 : nSubresultants_k1
            
            for i2 = 1 : 1 : nSubresultants_k2
                
                % Get the singular values
                vSingularValues = svd(arr_Sk1k2{i1,i2});
                arr_SingularValues{i1, i2} = vSingularValues;
                mat_MinimumSingularValues(i1,i2) = min(min(arr_SingularValues{i1, i2}));
                
            end
        end
        
        % Plot Graphs
        if(SETTINGS.PLOT_GRAPHS)
            
            %plotSingularValues(arr_SingularValues, limits_k1, limits_k2, limits_t1, limits_t2)
            plotMinimumSingularValues(mat_MinimumSingularValues, ...
                limits_k1, limits_k2, ...
                limits_t1, limits_t2)
            
        end
        
        % Set metric
        metric = mat_MinimumSingularValues;
        
        
        
    case 'Normalised Minimum Singular Values'
        
        
        
        % Initialise some vectors
        %arr_SingularValues = cell(nSubresultants_k1, nSubresultants_k2);
        arr_NormalisedSingularValues = cell(nSubresultants_k1, nSubresultants_k2);
        
        % Initialise matrix to store minimum singular values of SVD of each
        % S_{k1,k2}
        %mat_MinimumSingularValues = zeros(min(m1,n1)+1, min(m2,n2)+1);
        mat_NormalisedMinimumSingularValues = zeros(min(m1,n1)+1, min(m2,n2)+1);
        
        for i1 = 1:1:nSubresultants_k1
            
            for i2 = 1:1:nSubresultants_k2
                
                % Get the singular values
                
                vSingularValues = svd(arr_Sk1k2{i1,i2});
                
                arr_NormalisedSingularValues{i1, i2} = vSingularValues./ vSingularValues(1);
                
                mat_NormalisedMinimumSingularValues(i1, i2) = min(min(arr_NormalisedSingularValues{i1,i2}));
                
                
                
            end
        end
        
        % Plot Graphs
        if(SETTINGS.PLOT_GRAPHS)
            
            plotMinimumSingularValues(mat_NormalisedMinimumSingularValues, ...
                limits_k1, limits_k2, limits_t1, limits_t2)
            
        end
        
        % Set metric
        
        metric = mat_NormalisedMinimumSingularValues;
    case 'Residuals'
        
        error('Error : Code not yet developed for this branch')
        
    otherwise
        
        error('Error : Not a valid metric');
        
end


% 26/09/2017 - A more robust method is required here

% Compute the degree of the GCD
delta_x = diff(log10(metric),1,1);
delta_y = diff(log10(metric),1,2);


new_delta_y = [delta_y, zeros(min(m1,n1) + 1, 1)];
new_delta_x = [delta_x ; zeros(1, min(m2,n2) + 1)];

total_delta_x_y = new_delta_x + new_delta_y;
total_delta_x_y(1,1) = 0;

[M, I] = max(total_delta_x_y(:));
[r,c] = ind2sub(size(total_delta_x_y),I);

t1 = r - 1;
t2 = c - 1;

%vDelta_x = sum(delta_x,2);
%[~, index] = max(vDelta_x);
%t1 = index + lowerLimit_k1 -1;


%vDelta_y = sum(delta_y,1);
%[~, index] = max(vDelta_y);
%t2 = index + lowerLimit_k1 -1;


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

function [] = PlotAlphas(matrix_alpha, limits_k1, limits_k2, ...
    limits_t1, limits_t2)
% % Inputs
%
% matMinimumSingularValues : (Matrix) Stores the minimum singular values of
% each Sylvester subresultant matrix S_{k1,k2} for k1 =
% lowerlim,...,upperlim and k2 = lowerlim,...,upperlim
%
% limits_k1 : [Int Int]
%
% limits_k2 : [Int Int]
%
% limits_t1 : [Int Int]
%
% limits_t2 : [Int Int]
%
% rank_range : [Float Float]

global SETTINGS

% Get upper and lower limit
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);

%
lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

%
v_i1 = lowerLimit_k1 : 1 : upperLimit_k1;
v_i2 = lowerLimit_k2 : 1 : upperLimit_k2;




try
    [X,Y] = meshgrid(v_i1,v_i2);
    
    figure_name = sprintf('alpha : %s', SETTINGS.SYLVESTER_MATRIX_VARIANT);
    figure('name',figure_name)
    hold on
    %title('Minimum Singular Values');
    
    surf(X, Y, log10(matrix_alpha)');
    
    % Labels
    xlabel('$k_{1}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$k_{2}$', 'Interpreter', 'latex', 'FontSize', 20);
    zlabel('$\log_{10} \left( \alpha_{k_{1}, k_{2}} \right)$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Set view angle
    az = -30;
    el = 20;
    view(az, el);
    
    
    % Display
    grid on
    box on
    
    % Set location of window and size
    m_left = 100;
    m_bottom = 100;
    m_width = 600;
    m_height = 600;
    
    set(gcf, 'Position', [m_left, m_bottom, m_width, m_height]);
    
    % Position of figure within window
    myplot = gca;
    myval_side = 0.12;
    myval_base = 0.10;
    set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])
    
    
    hold off
catch
    
end







end


