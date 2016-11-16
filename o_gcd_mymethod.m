function [fxy_o, gxy_o, dxy_o, uxy_o, vxy_o, t, t1, t2] = ...
    o_gcd_mymethod(fxy, gxy, m, n, t_limits)
% o1(fxy_matrix,gxy_matrix,m,n)
%
% Given two input polynomials, calculate the GCD and its degree structure
%
% % Inputs.
%
%
% fxy_matrix : Matrix of coefficients of f(x,y) in the Bernstein basis.
%
% gxy_matrix : Matrix of coefficients of g(x,y) in the Bernstein basis
%
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomial g(x,y)
%
% t_limits : Limits on degree of GCD.
%
% % Outputs
%
% fxy : 
%
% gxy : 
%
% dxy : 
%
% uxy : 
%
% vxy :
% 
% t
% 
% t1 : 
% 
% t2 :




% Get Degree by first finding the total degree, then obtain t1 and t2

% Get total degreee of the GCD
[t, ~, ~] = GetGCDDegree_Total(fxy, gxy, m, n,t_limits);
fprintf([mfilename ' : ' sprintf('Degree of GCD : %i \n',t)])

% Get degree of the GCD with respect to x and y (t_{1} and t_{2})
[t1,t2,lambda,mu,alpha, th1,th2] = GetGCDDegree_Relative(fxy,gxy);
fprintf([mfilename ' : ' sprintf('Degree of GCD : t1 = %i, t2 = %i \n',t1,t2)])

fxy_n = fxy./lambda;
gxy_n = gxy./mu;

fww_n = GetWithThetas(fxy_n ,th1,th2);
gww_n = GetWithThetas(gxy_n ,th1,th2);
a_gww_n = alpha.* gww_n;

% Get optimal column
St1t2 = BuildDTQ(fww_n,a_gww_n,t1,t2);
idx_col = GetOptimalColumn(St1t2);

% % Get coefficients of u(x,y) and v(x,y).
[fxy_lr,gxy_lr,uxy_lr,vxy_lr,alpha_lr,th1_lr,th2_lr] = ...
    LowRankApproximation(fxy_n, gxy_n, alpha, th1, th2, m,n,t,t1, t2, idx_col);


% % Get coefficients of d(x,y)
[fxy_lra,gxy_lra,uxy_lra,vxy_lra,dxy_lra,alpha_lra,th1_lra,th2_lra] = ...
    APF(fxy_lr, gxy_lr, uxy_lr, vxy_lr, m, n, t, t1, t2, alpha_lr, th1_lr, th2_lr);


% % 
% Outputs
fxy_o = fxy_lra;
gxy_o = gxy_lra;
uxy_o = uxy_lra;
vxy_o = vxy_lra;
dxy_o = dxy_lra;







% % Compare Singular values of S(f(x,y),g(x,y)) and S(f+\delta f, g+ \delta g)
S_unproc = BuildDTQ(fxy,gxy,0,0);

% Preprocess f(x,y) to obtain f(w,w)
fw_preproc = GetWithThetas(fxy,th1,th2);

% Preprocess g(x,y) to obtain g(w,w)
gw_preproc = GetWithThetas(gxy,th1,th2);

% Build the 0th subresultants
S_preproc = BuildDTQ(fw_preproc,alpha.*gw_preproc,0,0);

% Build the Sylvester matrix with preprocessed coefficients
fw_lr = GetWithThetas(fxy,th1,th2);
gw_lr = GetWithThetas(gxy,th1,th2);
S_lowrank = BuildDTQ(fw_lr,alpha.*gw_lr,0,0);

% Get singular values of S(f,g), S(f,g) and S(f,g)
vSingularValues_unproc = svd(S_unproc);
vSingularValues_preproc = svd(S_preproc);
vSingularValues_lowrank = svd(S_lowrank);

global SETTINGS
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf([mfilename ' : ' 'Sylvester Matrices']);
        figure('name',figure_name)
        hold on
        xlabel('k')
        ylabel('log_{10}')
        plot(log10(vSingularValues_unproc),'-s','DisplayName','S(f(x,y),g(x,y))')
        plot(log10(vSingularValues_preproc),'-s','DisplayName','S(f(\omega,\omega),g(\omega,\omega))')
        plot(log10(vSingularValues_lowrank),'-s','DisplayName','Low Rank Approx')
        legend(gca,'show');
        hold off
    case 'n'
    otherwise
        error('err')
end
end




