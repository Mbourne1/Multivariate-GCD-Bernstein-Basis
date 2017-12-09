function [t1, t2, GM_fxy, GM_gxy, GM_hxy, alpha, beta, gamma, th1, th2, rank_range] = ...
    GetGCDDegree_Bivariate_3Polys(fxy, gxy, hxy, limits_t1, limits_t2, rank_range)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y)
% 
% % Inputs.
%
% fxy : (Matrix) Coefficient matrix of polynomial f(x,y)
%
% gxy : (Matrix) Coefficient matrix of polynomial g(x,y)
%
% hxy : (Matrix) Coefficient matrix of polynomial h(x,y)
%
% limits_t1 : (Int Int) :
%
% limits_t2 : (Int Int) : 
%
% % Outputs
%
% t1 : (Int) Degree of d(x,y) with respect to x
%
% t2 : (Int) Degree of d(x,y) with respect to y
%
% GM_fx : (Float) Geometric mean of entries in the first partition of the Sylvester
% subresultatn matrix S_{t1,t2}(f,g)
%
% GM_gx : (Float) Geometric mean of entries in the second partition of the
% Sylvester subresultant matrix S_{t1,t2}(f,g)
%
% alpha : (Float) Optimal value of alpha
%
% th1 : (Float) Optimal value of theta_{1}
%
% th2 : (Float) Optimal value of theta_{2}
%
% rank_range [(Float) (Float)]


%
global SETTINGS


% Get the degree structure of polynomial f(x,y), g(x,y) and h(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);

% Set my limits in the computation of the degree of the GCD.
limits_k1 = [0 min([m1,n1,o1])];
limits_k2 = [0 min([m2,n2,o2])];

%
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);

%
lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

% 
rank_range_low = rank_range(1);
rank_range_high = rank_range(2);


% Get number of Sylvester subresultant matrices to be considered.
nSubresultants_k1 = upperLimit_k1 - lowerLimit_k1 + 1;
nSubresultants_k2 = upperLimit_k2 - lowerLimit_k2 + 1;

% Initialise some vectors
arrSk1k2 = cell(nSubresultants_k1, nSubresultants_k2);
matAlpha = zeros(nSubresultants_k1, nSubresultants_k2);
matBeta = zeros(nSubresultants_k1, nSubresultants_k2);
matTheta1 = zeros(nSubresultants_k1, nSubresultants_k2);
matTheta2 = zeros(nSubresultants_k1, nSubresultants_k2);

matGM_fxy = zeros(nSubresultants_k1, nSubresultants_k2);
matGM_gxy = zeros(nSubresultants_k1, nSubresultants_k2);
matGM_hxy = zeros(nSubresultants_k1, nSubresultants_k2);

% For every Sylvester subresultant

for i1 = 1 : 1 : nSubresultants_k1
    
    for i2 = 1:1:nSubresultants_k2
        
        % Get k_{1} and k_{2}, the index of the S_{k1,k2} matrix
        k1 = lowerLimit_k1 + (i1 - 1); 
        k2 = lowerLimit_k2 + (i2 - 1);
        
        % Preprocess polynomials f(x,y), g(x,y) and h(x,y)
        [GM_fxy, GM_gxy, GM_hxy, alpha, beta, gamma, th1, th2] = ...
            Preprocess_Bivariate_3Polys(fxy, gxy, hxy, k1, k2);
        
        % Store the values from preprocessing
        matGM_fxy(i1,i2) = GM_fxy;
        matGM_gxy(i1,i2) = GM_gxy;
        matGM_hxy(i1,i2) = GM_hxy;
        matTheta1(i1,i2) = th1;
        matTheta2(i1,i2) = th2;
        matAlpha(i1,i2) = alpha;
        matBeta(i1,i2) = beta;
        
        % Divide f(x,y), g(x,y) and h(x,y) by geometric mean
        fxy_n = fxy ./ GM_fxy;
        gxy_n = gxy ./ GM_gxy;
        hxy_n = hxy ./ GM_hxy;
        
        
      
        
        
        
        % Get f(w1,w2), g(w1,w2) and h(w1,w2)
        lambda_fww = alpha .* GetWithThetas(fxy_n, th1, th2);
        mu_gww = beta .* GetWithThetas(gxy_n, th1, th2);
        rho_hww = gamma .* GetWithThetas(hxy_n, th1, th2);
        
        % Build the k1,k2 subresultant matrix
        arrSk1k2{i1,i2} = BuildSubresultant_Bivariate_3Polys(lambda_fww, mu_gww, rho_hww, k1, k2);
        
        if i1 == 1 && i2 == 1
        
            arr_Polys = {fxy, gxy, hxy, lambda_fww, mu_gww, rho_hww};
            arr_Names = {'$f(x,y)$', '$g(x,y)$', '$h(x,y)$', ...
                '$\lambda \tilde{f}(\omega_{1}, \omega_{2})$',...
                '$\mu \tilde{g}(\omega_{1}, \omega_{2})$',...
                '$\rho \tilde{h}(\omega_{1}, \omega_{2})$'};
            PlotCoefficients(arr_Polys, arr_Names);
            
        end
        
        
    end
    
end

% Metric used to compute the degree of the GCD
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals
switch SETTINGS.RANK_REVEALING_METRIC
    
    case 'R1 Row Norms'
        
        error('error : Code not complete')
        
    case 'R1 Row Diagonals'
        
        error('error : Code not complete')
        
    case 'Minimum Singular Values'
        
        matMinimumSingularValues = zeros(nSubresultants_k1, nSubresultants_k2);
        matNormalisedMinimumSingularValues = zeros(nSubresultants_k1, nSubresultants_k2);
        
        for i1 = 1 : 1 : nSubresultants_k1
            for i2 = 1 : 1 : nSubresultants_k2
                
                % Get the singular values of S_{k1,k2}
                vSingularValues = svd(arrSk1k2{i1,i2});
                
                vNormalisedSingularValues = vSingularValues./vSingularValues(1);
                
                % Get minimum singular value
                matMinimumSingularValues(i1,i2) = min(vSingularValues);
                matNormalisedMinimumSingularValues(i1,i2) = min(vNormalisedSingularValues);
                
            end
        end
        
        if (SETTINGS.PLOT_GRAPHS)
            
            %plotSingularValues(arrSingularValues, myLimits_t1, myLimits_t2, limits_t1, limits_t2);
            %plotMinimumSingularValues( matMinimumSingularValues, limits_k1, limits_k2, limits_t1, limits_t2, rank_range);
            
            plotMinimumSingularValues( matNormalisedMinimumSingularValues, limits_k1, limits_k2, limits_t1, limits_t2, rank_range);
            
        end
        
        vMetric = log10(matNormalisedMinimumSingularValues);
        
        
        
    case 'Residuals'
        
        error('error : Code not complete')
        
    otherwise
        
        error('Not a valid method')
        
end



% Compute the degree of the GCD
delta_x = diff(vMetric,1,1);
vDelta_x = sum(delta_x,2);
[~, index1] = max(vDelta_x);
t1 = index1 + lowerLimit_k1 -1;

delta_y = diff(vMetric,1,2);
vDelta_y = sum(delta_y,1);
[~, index2] = max(vDelta_y);
t2 = index2 + lowerLimit_k1 -1;


% Establish new rank range
x1 = index1;
x2 = index2;
rank_range_low = vMetric(x1, x2);
rank_range_high = vMetric(x1 + 1, x2 + 1);

rank_range = [rank_range_low rank_range_high];

% Outputs
GM_fxy = matGM_fxy(t1 - lowerLimit_k1 +1, t2 - lowerLimit_k2 +1);
GM_gxy = matGM_gxy(t1 - lowerLimit_k1 +1, t2 - lowerLimit_k2 +1);
GM_hxy = matGM_hxy(t1 - lowerLimit_k1 +1, t2 - lowerLimit_k2 +1);
alpha = matAlpha(t1 - lowerLimit_k1 +1, t2 - lowerLimit_k2 +1);
beta = matBeta(t1 - lowerLimit_k1 +1, t2 - lowerLimit_k2 +1);
th1 = matTheta1(t1 - lowerLimit_k1 +1, t2 - lowerLimit_k2 +1);
th2 = matTheta2(t1 - lowerLimit_k1 +1, t2 - lowerLimit_k2 +1);


LineBreakMedium()
fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n', t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n', t2)])
LineBreakMedium()

end


function PlotCoefficients(arr_Polys, arr_Names)

nPolys = length(arr_Polys);

figure()

hold on

for i = 1 : 1 : nPolys
   
    fxy = arr_Polys{i};
    fxy = GetAsVector(fxy);
    vec_x = 1 : 1: length(fxy);
    fx_str = arr_Names{i};
    plot(vec_x, log10(fxy),'DisplayName',fx_str);
        
end

l = legend(gca,'show')
set(l, 'Interpreter','latex');
set(l, 'FontSize',20);
box on
grid on


end



