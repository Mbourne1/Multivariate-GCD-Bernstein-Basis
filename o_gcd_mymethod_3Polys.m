function [fxy_o, gxy_o, hxy_o, dxy_o, uxy_o, vxy_o, wxy_o, t, t1, t2] = ...
    o_gcd_mymethod_3Polys(fxy, gxy, hxy, m, n, o, t_limits)
% o1(fxy_matrix,gxy_matrix,m,n)
%
% Given two input polynomials, calculate the GCD and its degree structure
%
% % Inputs.
%
%
% [fxy, gxy, hxy] : Matrix of coefficients of f(x,y) g(x,y) and h(x,y) in 
% the Bernstein basis.
%
% [m, n, o] : Total degree of polynomial f(x,y), g(x,y) and h(x,y)
%
% t_limits : Limits on degree of GCD.
%
% % Outputs
%
% [fxy, gxy, hxy] : Coefficients of the polynomial f(x,y)
%
% dxy : Coefficients of the polynomial d(x,y)
%
% [uxy, vxy, wxy] : Coefficients of the polynomial u(x,y)
%
% t : Total degree of the GCD
% 
% t1 : Degree of the GCD with respect to x
% 
% t2 : Degree of the GCD with respect to y

% Get Degree by first finding the total degree, then obtain t1 and t2

% Get total degreee of the GCD
[t, ~, ~] = GetGCDDegree_Total_3Polys(fxy, gxy, hxy, m, n, o, t_limits);

fprintf([mfilename ' : ' sprintf('Degree of GCD : %i \n',t)])

% Get degree of the GCD with respect to x and y (t_{1} and t_{2})
[t1, t2, GM_fx, GM_gx, alpha, th1, th2] = GetGCDDegree_Relative_3Polys(fxy, gxy, hxy);

fprintf([mfilename ' : ' sprintf('Degree of GCD : t1 = %i, t2 = %i \n',t1,t2)])

% Normalise f(x,y) and g(x,y)
fxy_n = fxy./GM_fx;
gxy_n = gxy./GM_gx;
hxy_n = hxy./GM_hx;

% Preprocess f(x,y) and g(x,y)
fww_n = GetWithThetas(fxy_n, th1, th2);
gww_n = GetWithThetas(gxy_n, th1, th2);
hww_n = GetWithThetas(hxy_n, th1, th2);

a_gww_n = alpha.* gww_n;

% Get optimal column
St1t2 = BuildDTQ_3Polys(fww_n,a_gww_n, hww_n, t1, t2);
idx_col = GetOptimalColumn(St1t2);

% % Get coefficients of u(x,y) and v(x,y).
[fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    LowRankApproximation(fxy_n, gxy_n, hxy_n, alpha, th1, th2, m, n, o, t, t1, t2, idx_col);


% % Get coefficients of d(x,y)
[fxy_lra, gxy_lra, uxy_lra, vxy_lra, dxy_lra, alpha_lra, th1_lra, th2_lra] = ...
    APF(fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, m, n, o, t, t1, t2, alpha_lr, th1_lr, th2_lr);


% % 
% Outputs
fxy_o = fxy_lra;
gxy_o = gxy_lra;
hxy_o = hxy_lra;

uxy_o = uxy_lra;
vxy_o = vxy_lra;
wxy_o = wxy_lra;

dxy_o = dxy_lra;







end




