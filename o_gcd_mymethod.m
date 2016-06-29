function [fxy,gxy,dxy,uxy, vxy,t,t1,t2] = o_gcd_mymethod(fxy,gxy,...
    m,n,t_limits)
% o1(fxy_matrix,gxy_matrix,m,n)
%
% Given two input polynomials, calculate the GCD and its degree structure
%
% % Inputs.
%
%
% fxy_matrix : Matrix of coefficients of f(x,y)
%
% gxy_matrix : Matrix of coefficients of g(x,y)
%
% m : Total degree of polynomial f(x,y).
%
% n : Total degree of polynomial g(x,y).

% %
% % Get the degree of the GCD d(x,y)
degree_calc_method = 'respective';

input_fxy = fxy;
input_gxy = gxy;

% Get Degree by first finding the total degree, then obtain t1 and t2

% Get total degreee
%[t_old, ~, ~] = GetGCDDegree_Total(fxy_matrix, gxy_matrix,m,n);
[t_new, ~, ~] = GetGCDDegree_Total2(fxy, gxy,m,n, t_limits);


t = t_new;
fprintf([mfilename ' : ' sprintf('Degree of GCD : %i \n',t)])


% Get degree t1 and t2
[t1,t2,lambda,mu,alpha, th1,th2] = GetGCDDegree_Relative(fxy,gxy,m,n,t);

fprintf([mfilename ' : ' sprintf('Degree of GCD : t1 = %i, t2 = %i \n',t1,t2)])

fxy_matrix_n = fxy./lambda;
gxy_matrix_n = gxy./mu;

fww_matrix_n = GetWithThetas(fxy_matrix_n ,th1,th2);
gww_matrix_n = GetWithThetas(gxy_matrix_n ,th1,th2);
a_gww_matrix_n = alpha.* gww_matrix_n;

% Get Optimal column for removal from S_{t_{1},t_{2}}
opt_col = GetOptimalColumn(fww_matrix_n,a_gww_matrix_n,t1,t2);


%
[fxy,gxy,alpha,th1,th2] = LowRankApproximation...
    (fxy_matrix_n,gxy_matrix_n,alpha,th1,th2,t1,t2,opt_col);



% % Get Quotients u(x,y) and v(x,y)
% Calc method is either total or respective



switch degree_calc_method
    case 'respective'
        
        fww = GetWithThetas(fxy_matrix_n,th1,th2);
        gww = GetWithThetas(gxy_matrix_n,th1,th2);
        a_gww = alpha.* gww;
        
        [uww, vww] ...
            = GetQuotients(fww, a_gww,t1,t2);
        
        uxy = GetWithoutThetas(uww,th1,th2);
        vxy = GetWithoutThetas(vww,th1,th2);
        
    case 'total'
        [uxy, vxy] ...
            = GetQuotients_total(fww_matrix_n, alpha.*gww_matrix_n,m,n,t);
    otherwise
        
end




% % Get the GCD
% % Get d(x,y) from the polynomials u(x,y) and v(x,y).
switch degree_calc_method
    case 'respective'
        dww_calc_matrix = GetGCD_Coefficients(uww,vww,...
            fww,alpha.*gww,t1,t2);
    case 'total'
        dww_calc_matrix = GetGCD_Coefficients_total(uww,vww,...
            fww_matrix,gww_matrix,m,n,t);
    otherwise
        error('error')
end

dxy = GetWithoutThetas(dww_calc_matrix,th1,th2);

% % Compare Singular values of S(f(x,y),g(x,y)) and S(f+\delta f, g+ \delta g)
S_unproc = BuildDTQ(input_fxy,input_gxy,0,0);

% Preprocess f(x,y) to obtain f(w,w)
fw_preproc = GetWithThetas(input_fxy,th1,th2);

% Preprocess g(x,y) to obtain g(w,w)
gw_preproc = GetWithThetas(input_gxy,th1,th2);

% Build the 0th subresultants
S_preproc = BuildDTQ(fw_preproc,alpha.*gw_preproc,0,0);

fw_lr = GetWithThetas(fxy,th1,th2);
gw_lr = GetWithThetas(gxy,th1,th2);
S_lowrank = BuildDTQ(fw_lr,alpha.*gw_lr,0,0);

vSingularValues_unproc = svd(S_unproc);
vSingularValues_preproc = svd(S_preproc);
vSingularValues_lowrank = svd(S_lowrank);

figure_name = sprintf('%s - Sylvester Matrices',mfilename);
figure('name',figure_name)
hold on
plot(log10(vSingularValues_unproc),'DisplayName','S(f(x,y),g(x,y))')
plot(log10(vSingularValues_preproc),'DisplayName','S(f(\omega,\omega),g(\omega,\omega))')
plot(log10(vSingularValues_lowrank),'DisplayName','Low Rank Approx')
legend(gca,'show');
hold off

end




