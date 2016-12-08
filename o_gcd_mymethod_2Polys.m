function [fxy_o, gxy_o, dxy_o, uxy_o, vxy_o, t, t1, t2] = ...
    o_gcd_mymethod_2Polys(fxy, gxy, m, n, t_limits)
% o1(fxy_matrix,gxy_matrix,m,n)
%
% Given two input polynomials, calculate the GCD and its degree structure
%
% % Inputs.
%
%
% [fxy, gxy] : Matrix of coefficients of f(x,y) and g(x,y) in the Bernstein basis.
%
% [m, n] : Total degree of polynomial f(x,y)
%
% t_limits : Limits on degree of GCD.
%
% % Outputs
%
% [fxy_o, gxy_o] : Coefficients of the polynomial f(x,y) and g(x,y)
%
% dxy_o : Coefficients of the polynomial d(x,y)
%
% [uxy_o, vxy_o] : Coefficients of the polynomial u(x,y) and v(x,y)
%
% t : Total degree of the GCD
% 
% [t1, t2] : Degree of the GCD with respect to x and y
% 


% Get Degree by first finding the total degree, then obtain t1 and t2

% Get total degreee of the GCD
[t, ~, ~] = GetGCDDegree_Total_2Polys(fxy, gxy, m, n,t_limits);
fprintf([mfilename ' : ' sprintf('Degree of GCD : %i \n',t)])

% Get degree of the GCD with respect to x and y (t_{1} and t_{2})
[t1,t2,lambda,mu,alpha, th1,th2] = GetGCDDegree_Relative_2Polys(fxy, gxy);
fprintf([mfilename ' : ' sprintf('Degree of GCD : t1 = %i, t2 = %i \n',t1,t2)])

% Normalise f(x,y) and g(x,y)
fxy_n = fxy./lambda;
gxy_n = gxy./mu;

% Preprocess f(x,y) and g(x,y)
fww_n = GetWithThetas(fxy_n ,th1,th2);
gww_n = GetWithThetas(gxy_n ,th1,th2);
a_gww_n = alpha.* gww_n;

% Get optimal column
St1t2 = BuildDTQ_2Polys(fww_n,a_gww_n,t1,t2);
idx_col = GetOptimalColumn(St1t2);

% % Get coefficients of u(x,y) and v(x,y).
[fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    LowRankApproximation_2Polys(fxy_n, gxy_n, alpha, th1, th2, m, n, t, t1, t2, idx_col);


% % Get coefficients of d(x,y)
[fxy_lra,gxy_lra,uxy_lra,vxy_lra,dxy_lra,alpha_lra,th1_lra,th2_lra] = ...
    APF_2Polys(fxy_lr, gxy_lr, uxy_lr, vxy_lr, m, n, t, t1, t2, alpha_lr, th1_lr, th2_lr);


% % 
% Outputs
fxy_o = fxy_lra;
gxy_o = gxy_lra;
uxy_o = uxy_lra;
vxy_o = vxy_lra;
dxy_o = dxy_lra;







end




