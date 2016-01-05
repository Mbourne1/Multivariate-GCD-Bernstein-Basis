function opt_col = getOptimalColumn(fxy_matrix,gxy_matrix,t1,t2,lambda,mu,opt_alpha,opt_theta_1,opt_theta_2)
%% Find Optimal column for removal from S_{t_{1},t_{2}}
% Given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized

% Inputs

% fxy_matrix : Matrix of coefficients of polynomial f 
% gxy_matrix : Matrix of coefficients of polynomial g
% t1         : Degree of GCD with respect to x
% t2         : Degree of GCD with respect to y
% lambda     : Geometric mean of polynomial f
% mu         : Geometric mean of polynomial g
% opt_alpha  : Optimal value of alpha
% opt_theta_1: Optimal value of theta1
% opt_theta_2: Optimal value of theta2


%%

% Normalise f by geometric mean
fxy_matrix_n = fxy_matrix./lambda;

% Normalise g by geometric mean
gxy_matrix_n = gxy_matrix./mu;

% Get the degree of polynomial f(x,y)
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get the degree of polynomial g(x,y)
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

% Build the (t_{1},t_{2})-th subresultant
St1t2 = BuildSubresultant(fxy_matrix_n,gxy_matrix_n,t1,t2,opt_alpha, opt_theta_1, opt_theta_2);

%% Find Optimal column for removal from St
% given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized

[~,cols] = size(St1t2);

% QR Decomposition of the Sylvester Matrix S_{k}
[Qk,Rk] = qr(St1t2);
n = n1+n2;
for j=1:1:cols
    
    ck = St1t2(:,j);
    [Q,~] = qrdelete(Qk,Rk,j);
    cd = Q'*ck;
    d = cd(n+1:end,:);
    residuals_QR(j) = norm(d);
    
end

%Obtain the column for which the residual is minimal.
[~,opt_col] = min(log10(residuals_QR));
fprintf('Optimal column for removal is given by %i \n',opt_col)

end