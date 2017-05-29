function [fxy_o, gxy_o, hxy_o, dxy_o, uxy_o, vxy_o, wxy_o, t1, t2, rank_range] = ...
    o_gcd_mymethod_Bivariate_3Polys(fxy, gxy, hxy, limits_t1, limits_t2, rank_range)
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

%
[t1, t2, GM_fx, GM_gx, GM_hx, alpha, beta, th1, th2, rank_range] = GetGCDDegree_Bivariate_3Polys(fxy, gxy, hxy, limits_t1, limits_t2, rank_range);

fprintf([mfilename ' : ' sprintf('Degree of GCD : t1 = %i, t2 = %i \n', t1, t2)])

% Normalise f(x,y) and g(x,y)
fxy_n = fxy ./ GM_fx;
gxy_n = gxy ./ GM_gx;
hxy_n = hxy ./ GM_hx;

% Preprocess f(x,y) and g(x,y)
fww_n = GetWithThetas(fxy_n, th1, th2);
gww_n = GetWithThetas(gxy_n, th1, th2);
hww_n = GetWithThetas(hxy_n, th1, th2);

a_gww_n = alpha .* gww_n;
b_hww_n = beta .* hww_n;

% Get optimal column
St1t2 = BuildDTQ_Bivariate_3Polys(fww_n, a_gww_n, b_hww_n, t1, t2);

idx_col = GetOptimalColumn(St1t2);

% % Get coefficients of u(x,y) and v(x,y).
[fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, alpha_lr, beta_lr, th1_lr, th2_lr] = ...
    LowRankApproximation_Bivariate_3Polys(fxy_n, gxy_n, hxy_n, alpha, beta, th1, th2, t1, t2, idx_col);

% Get coefficients of d(x,y)
[fxy_lra, gxy_lra, hxy_lra, uxy_lra, vxy_lra, wxy_lra, dxy_lra, alpha_lra, th1_lra, th2_lra] = ...
    APF_Bivariate_3Polys(fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, t1, t2, alpha_lr, beta_lr, th1_lr, th2_lr);

% Outputs
fxy_o = fxy_lra;
gxy_o = gxy_lra;
hxy_o = hxy_lra;

uxy_o = uxy_lra;
vxy_o = vxy_lra;
wxy_o = wxy_lra;

dxy_o = dxy_lra;

end




