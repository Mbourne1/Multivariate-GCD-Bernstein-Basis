function [fxy_o, gxy_o, dxy_o, uxy_o, vxy_o, t, t1, t2] = ...
    o_gcd_mymethod_Bivariate_2Polys(fxy, gxy, m, n, t_limits)
% o1(fxy_matrix,gxy_matrix,m,n)
%
% Given two input polynomials, calculate the GCD and its degree structure
%
% % Inputs.
%
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) in the Bernstein basis.
%
% gxy : (Matrix) Coefficients of polynomial g(x,y) in the Bernstein basis.
%
% [m, n] : (Int) Total degree of polynomials f(x,y) and g(x,y)
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




%[t, ~, ~] = GetGCDDegree_Total_Bivariate_2Polys(fxy, gxy, m, n, t_limits);
t = 0;
fprintf([mfilename ' : ' sprintf('Degree of GCD : %i \n',t)])

% Get degree of the GCD with respect to x and y (t_{1} and t_{2})
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);

% Get upper and lower limits of degree of GCD
limits_t1 = [0 min(m1, n1)];
limits_t2 = [0 min(m2, n2)];


degree_method = 'All Subresultants';

switch degree_method
    case 'All Subresultants'
        
        [t1, t2, GM_fx, GM_gx, alpha, th1, th2] = ...
            GetGCDDegree_Bivariate_2Polys(fxy, gxy, limits_t1, limits_t2);
        
    case 'Linear Method'
        
        [t1, t2, GM_fx, GM_gx, alpha, th1, th2] = ...
            GetGCDDegree_Bivariate_2Polys_NewMethod(fxy, gxy, limits_t1, limits_t2);
        
end


fprintf([mfilename ' : ' sprintf('Degree of GCD : t1 = %i, t2 = %i \n',t1,t2)])

% Normalise f(x,y) and g(x,y)
fxy_n = fxy./ GM_fx;
gxy_n = gxy./ GM_gx;

% Preprocess f(x,y) and g(x,y)
fww_n = GetWithThetas(fxy_n, th1, th2);
gww_n = GetWithThetas(gxy_n, th1, th2);
a_gww_n = alpha.* gww_n;

% Get optimal column

St1t2 = BuildDTQ_Bivariate_2Polys(fww_n, a_gww_n, t1, t2);
idx_col = GetOptimalColumn(St1t2);

% % Get coefficients of u(x,y) and v(x,y).
[fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    LowRankApproximation_Bivariate_2Polys(fxy_n, gxy_n, alpha, th1, th2, m, n, t, t1, t2, idx_col);


% % Get coefficients of d(x,y)
[fxy_lra,gxy_lra,uxy_lra,vxy_lra,dxy_lra,alpha_lra,th1_lra,th2_lra] = ...
    APF_Bivariate_2Polys(fxy_lr, gxy_lr, uxy_lr, vxy_lr, m, n, t, t1, t2, alpha_lr, th1_lr, th2_lr);


% %
% Outputs
fxy_o = fxy_lra;
gxy_o = gxy_lra;
uxy_o = uxy_lra;
vxy_o = vxy_lra;
dxy_o = dxy_lra;







end




