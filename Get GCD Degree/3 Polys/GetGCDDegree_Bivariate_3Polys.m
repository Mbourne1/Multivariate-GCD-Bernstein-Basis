function [t1, t2, GM_fxy, GM_gxy, GM_hxy, alpha, th1, th2] = GetGCDDegree_Bivariate_3Polys(fxy, gxy, hxy, limits_t1, limits_t2)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y)
%
% % Inputs.
%
% [fxy, gxy, hxy] : Coefficient matrix of polynomial f(x,y) and g(x,y)
%
% % Outputs 
%
% t1 : Degree of d(x,y) with respect to x
% 
% t2 : Degree of d(x,y) with respect to y
%
% [GM_fx, GM_gx] : Geometric mean of entries in each partition of the
% Sylvester subresultant matrix S_{t1,t2}(f,g)
%
% alpha : Optimal value of alpha
%
% th1 : Optimal value of theta_{1}
%
% th2 : Optimal value of theta_{2}


%
global SETTINGS


% Get the degree structure of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get the degree structure of polynomial g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Get the degree structure of polynomial g(x,y)
[o1, o2] = GetDegree_Bivariate(hxy);

my_limits_t1 = [0 min([m1,n1,o1])];
my_limits_t2 = [0 min([m2,n2,o2])];

lowerLimit_t1 = my_limits_t1(1);
upperLimit_t1 = my_limits_t1(2);

lowerLimit_t2 = my_limits_t2(1);
upperLimit_t2 = my_limits_t2(2);

nSubresultants_k1 = upperLimit_t1 - lowerLimit_t1 + 1;
nSubresultants_k2 = upperLimit_t2 - lowerLimit_t2 + 1;


% Initialise some vectors
arrSingularValues = cell(nSubresultants_k1, nSubresultants_k2);

matAlpha = zeros(nSubresultants_k1, nSubresultants_k2);
matTheta1 = zeros(nSubresultants_k1, nSubresultants_k2);
matTheta2 = zeros(nSubresultants_k1, nSubresultants_k2);
matGM_fxy = zeros(nSubresultants_k1, nSubresultants_k2);
matGM_gxy = zeros(nSubresultants_k1, nSubresultants_k2);
matGM_hxy = zeros(nSubresultants_k1, nSubresultants_k2);

% For every row in the matrix
for i1 = 1:1:nSubresultants_k1
    
    for i2 = 1:1:nSubresultants_k2
    
    k1 = lowerLimit_t1 + (i1 - 1);
    k2 = lowerLimit_t2 + (i2 - 1);
    
    
    % Preprocessing
    %[GM_fx, GM_gx, alpha, th1, th2] = Preprocess(fxy, gxy, v_k1, v_k2);
    
    [GM_fxy, GM_gxy, GM_hxy, alpha, th1, th2] = Preprocess_Bivariate_3Polys(fxy, gxy, hxy, k1, k2);
    
    matGM_fxy(i1,i2) = GM_fxy;
    matGM_gxy(i1,i2) = GM_gxy;
    matGM_hxy(i1,i2) = GM_hxy;
    matTheta1(i1,i2) = th1;
    matTheta2(i1,i2) = th2;
    matAlpha(i1,i2) = alpha;
    
    
    % Divide f(x,y), g(x,y) and h(x,y) by geometric mean
    fxy_n = fxy ./ GM_fxy;
    gxy_n = gxy ./ GM_gxy;
    hxy_n = hxy ./ GM_hxy;
    
    % Get f(w1,w2), g(w1,w2) and h(w1,w2)
    fww = GetWithThetas(fxy_n, th1, th2);
    gww = GetWithThetas(gxy_n, th1, th2);
    hww = GetWithThetas(hxy_n, th1, th2);
    
    % Build the k1,k2 subresultant matrix
    Sk1k2 = BuildDTQ_Bivariate_3Polys(fww,gww,hww,k1,k2);
    
    % Get the singular values of S_{k1,k2}
    arrSingularValues{i1,i2} = svd(Sk1k2);
    
    
    end
    
end



% Metric used to compute the degree of the GCD
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals
switch SETTINGS.RANK_REVEALING_METRIC

    case 'R1 Row Norms'
        
    case 'R1 Row Diagonals'
        
    case 'Singular Values'
        
        matMinimumSingularValues = zeros(nSubresultants_k1, nSubresultants_k2);
        
        for i1 = 1:1:nSubresultants_k1
            for i2 = 1:1:nSubresultants_k2
           
                matMinimumSingularValues(i1,i2) = min(arrSingularValues{i1,i2});
                
            end
        end
        
        plotSingularValues(arrSingularValues, my_limits_t1, my_limits_t2);
        plotMinimumSingularValues(matMinimumSingularValues, my_limits_t1, my_limits_t2);
        
        
        metric = matMinimumSingularValues;
        
    case 'Residuals'
    
end



% Compute the degree of the GCD
delta_x = diff(log10(metric),1,1);
vDelta_x = sum(delta_x,2);
[~, index] = max(vDelta_x);
t1 = index + lowerLimit_t1 -1;

delta_y = diff(log10(metric),1,2);
vDelta_y = sum(delta_y,1);
[~, index] = max(vDelta_y);
t2 = index + lowerLimit_t1 -1;


% Outputs
GM_fxy = matGM_fxy(t1 - lowerLimit_t1 +1, t2 - lowerLimit_t2 +1);
GM_gxy = matGM_gxy(t1 - lowerLimit_t1 +1, t2 - lowerLimit_t2 +1);
GM_hxy = matGM_hxy(t1 - lowerLimit_t1 +1, t2 - lowerLimit_t2 +1);
alpha = matAlpha(t1 - lowerLimit_t1 +1, t2 - lowerLimit_t2 +1);
th1 = matTheta1(t1 - lowerLimit_t1 +1, t2 - lowerLimit_t2 +1);
th2 = matTheta2(t1 - lowerLimit_t1 +1, t2 - lowerLimit_t2 +1);

LineBreakMedium()
fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n',t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n',t2)])
LineBreakMedium()

end




