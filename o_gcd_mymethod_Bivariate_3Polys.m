function [fxy_o, gxy_o, hxy_o, dxy_o, uxy_o, vxy_o, wxy_o, t1, t2, rank_range] = ...
    o_gcd_mymethod_Bivariate_3Polys(fxy, gxy, hxy, limits_t1, limits_t2, rank_range, degree_method)
% o_gcd_mymethod_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, t_limits)
%
% Given two input polynomials, calculate the GCD and its degree structure
%
% % Inputs.
%
%
% fxy : (Matrix) Coefficients of f(x,y)
%
% gxy : (Matrix) Coefficients of g(x,y)
%
% hxy : (Matrix) Coefficients of h(x,y)
%
% limits_t1 : [(Int) (Int)] Limits on degree of GCD with respect to x
% 
% limits_t2 : [(Int) (Int)] Limits on degree of GCD with respect to y
%
% rank_range : [(Float) (Float)]
%
% % Outputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of the polynomial h(x,y)
%
% dxy : (Matrix) Coefficients of the polynomial d(x,y), the greatest common
% divisor of f(x,y), g(x,y) and h(x,y).
%
% uxy : (Matrix) Coefficients of the polynomial u(x,y), the quotient
% polynomial such that f(x,y)/u(x,y) = d(x,y)
%
% vxy : (Matrix) Coefficients of the polynomial u(x,y), the quotient
% polynomial such that g(x,y)/v(x,y) = d(x,y)
%
% wxy : : (Matrix) Coefficients of the polynomial u(x,y), the quotient
% polynomial such that h(x,y)/w(x,y) = d(x,y)
%
% t : (Int) Total degree of the GCD
% 
% t1 :(Int) Degree of the GCD with respect to x
% 
% t2 :(Int) Degree of the GCD with respect to y


switch degree_method
    case 'All Subresultants'
        
        [t1, t2, GM_fxy, GM_gxy, GM_hxy, lambda, mu, rho, th1, th2, rank_range] ...
            = GetGCDDegree_Bivariate_3Polys(fxy, gxy, hxy, limits_t1, limits_t2, rank_range);

        
    case 'Linear Method'
        
        %[t1, t2, GM_fxy, GM_gxy, alpha, th1, th2] = ...
        %    GetGCDDegree_Bivariate_2Polys_Linear_Naive(fxy, gxy, limits_t1, limits_t2);
        
    case 'Total'    
        
         
        [t1_possible, t2_possible, lambda, mu, rho,  th1, th2] = GetGCDDegree_Bivariate_3Polys_WithDegreeElevation(fxy, gxy, hxy);
        
       
            
        candidate_t1_stage_a = t1_possible;
        [candidate_t1_stage_a, candidate_t2_stage_a, GM_fxy, GM_gxy, GM_hxy, lambda, mu, rho, th1, th2] = ...
            GetGCDDegree_Bivariate_3Polys_Linear(fxy, gxy, hxy, candidate_t1_stage_a, 'x');
        
        
        candidate_t2_stage_b = t2_possible;
        [candidate_t1_stage_b, candidate_t2_stage_b, GM_fxy, GM_gxy, GM_hxy, lambda, mu, rho, th1, th2] = ...
            GetGCDDegree_Bivariate_3Polys_Linear(fxy, gxy, hxy, candidate_t2_stage_b, 'y');
        
        fprintf('First Pair a: %i %i \n', candidate_t1_stage_a, candidate_t2_stage_a);
        fprintf('First Pair b: %i %i \n', candidate_t1_stage_b, candidate_t2_stage_b);
        
        
        t1 = candidate_t1_stage_a;
        t2 = candidate_t2_stage_a;
        
    otherwise 
        error('Not valid method. Must be *All Subresultants *Linear Method *Total')
        
end
%

fprintf([mfilename ' : ' sprintf('Degree of GCD : t1 = %i, t2 = %i \n', t1, t2)])

% Normalise f(x,y) and g(x,y)
fxy_n = fxy ./ GM_fxy;
gxy_n = gxy ./ GM_gxy;
hxy_n = hxy ./ GM_hxy;

% Preprocess f(x,y) and g(x,y)
a_fww = lambda .* GetWithThetas(fxy_n, th1, th2);
b_gww = mu .* GetWithThetas(gxy_n, th1, th2);
gamma_hww = rho.* GetWithThetas(hxy_n, th1, th2);


% Get optimal column
St1t2 = BuildSubresultant_Bivariate_3Polys(a_fww, b_gww, gamma_hww, t1, t2);
idx_col = GetOptimalColumn(St1t2);

% % Get coefficients of u(x,y) and v(x,y).
[fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, alpha_lr, beta_lr, gamma_lr, th1_lr, th2_lr] = ...
    LowRankApproximation_Bivariate_3Polys(fxy_n, gxy_n, hxy_n, lambda, mu, rho, th1, th2, t1, t2, idx_col);

% Get coefficients of d(x,y)
[fxy_lra, gxy_lra, hxy_lra, uxy_lra, vxy_lra, wxy_lra, dxy_lra, lambda_lra, mu_lra, rho_lra, th1_lra, th2_lra] = ...
    APF_Bivariate_3Polys(fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, t1, t2, alpha_lr, beta_lr, gamma_lr, th1_lr, th2_lr);

% Outputs
fxy_o = fxy_lra;
gxy_o = gxy_lra;
hxy_o = hxy_lra;

uxy_o = uxy_lra;
vxy_o = vxy_lra;
wxy_o = wxy_lra;

dxy_o = dxy_lra;

end




