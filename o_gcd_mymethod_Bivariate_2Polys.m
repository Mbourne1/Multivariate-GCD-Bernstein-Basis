function [fxy_o, gxy_o, dxy_o, uxy_o, vxy_o, t1, t2] = ...
    o_gcd_mymethod_Bivariate_2Polys(fxy, gxy, limits_t1, limits_t2, degree_method)
% o1(fxy_matrix,gxy_matrix,m,n)
%
% Given two input polynomials, calculate the GCD and its degree structure
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) in the Bernstein basis.
%
% gxy : (Matrix) Coefficients of polynomial g(x,y) in the Bernstein basis.
%
% m : (Int) Total degree of polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)
%
% limits_t1 : [(Int) (Int)] Limits on degree of GCD.
%
% limits_t2 :
%
% % Outputs
%
% fxy_o : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy_o : (Matrix) Coefficients of the polynomial g(x,y)
%
% dxy_o : (Matrix) Coefficients of the polynomial d(x,y)
%
% uxy_o : (Matrix) Coefficients of the polynomial u(x,y)
%
% vxy_o : (Matrix) Coefficients of the polynomial v(x,y)
%
% t : (Int) Total degree of the GCD
%
% t1 : (Int) Degree of GCD with respect to x
%
% t2 : (Int) Degree of GCD with respect to y


% Get Degree by first finding the total degree, then obtain t1 and t2

% 'All Subresultants Method' : Create all subresultants, get rank revealing
% metric for each S_{k1,k2} then compute the GCD by taking the max
% difference in both x and y.

% 'Linear Method' : Fix k_{1} = 1, produce all subresultants S_{1,k2}, compute
% t2, then fix k2 = t2 and compute all subresultants S_{k1,t2}


switch degree_method
    case 'All Subresultants'
        
        [t1, t2, GM_fx, GM_gx, alpha, th1, th2] = ...
            GetGCDDegree_Bivariate_2Polys(fxy, gxy, limits_t1, limits_t2);
        
    case 'Linear Method' % Dont use this method
        
        [t1, t2, GM_fx, GM_gx, alpha, th1, th2] = ...
            GetGCDDegree_Bivariate_2Polys_Linear(fxy, gxy, limits_t1, limits_t2);
        
    case 'Total'
        
      
        [t1_candidate, t2_candidate, alpha, th1, th2] = ... 
            GetGCDDegree_Bivariate_2Polys_WithDegreeElevation(fxy, gxy);
        
        
        % Assume t_{1} candidate is correct and compute t_{2}
        [t1_a, t2_a, GM_fx_a, GM_gx_a, alpha_a, th1_a, th2_a] = ...
            GetGCDDegree_Bivariate_2Polys_Linear(fxy, gxy, t1_candidate, 'x');
        
        % Assume t_{2} candidate is correct and compute t_{1}
        [t1_b, t2_b, GM_fx_b, GM_gx_b, alpha_b, th1_b, th2_b] = ...
            GetGCDDegree_Bivariate_2Polys_Linear(fxy, gxy, t2_candidate, 'y');
        
        % Note this is not correct
        t1 = t1_a;
        t2 = t2_a;
        GM_fx = GM_fx_a;
        GM_gx = GM_gx_a;
        alpha = alpha_a;
        th1 = th1_a;
        th2 = th2_a;
        
        
    otherwise
        error('Not valid method - Method must be *All Subresultants, *Linear Method or *Total')
        
end


fprintf([mfilename ' : ' sprintf('Degree of GCD : t1 = %i, t2 = %i \n',t1,t2)])

% Normalise f(x,y) and g(x,y)
fxy_n = fxy ./ GM_fx;
gxy_n = gxy ./ GM_gx;

% Preprocess f(x,y) and g(x,y)
fww_n = GetWithThetas(fxy_n, th1, th2);
gww_n = GetWithThetas(gxy_n, th1, th2);
a_gww_n = alpha .* gww_n;

% Get optimal column for removal
St1t2 = BuildSubresultant_Bivariate_2Polys(fww_n, a_gww_n, t1, t2);
idx_col = GetOptimalColumn(St1t2);

% % Get coefficients of u(x,y) and v(x,y).
[fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    LowRankApproximation_Bivariate_2Polys(fxy_n, gxy_n, alpha, th1, th2, t1, t2, idx_col);


% % Get coefficients of d(x,y)
[fxy_lra, gxy_lra, uxy_lra, vxy_lra, dxy_lra, alpha_lra, th1_lra, th2_lra] = ...
    APF_Bivariate_2Polys(fxy_lr, gxy_lr, uxy_lr, vxy_lr, t1, t2, alpha_lr, th1_lr, th2_lr);


% %
% Outputs
fxy_o = fxy_lra;
gxy_o = gxy_lra;
uxy_o = uxy_lra;
vxy_o = vxy_lra;
dxy_o = dxy_lra;







end




