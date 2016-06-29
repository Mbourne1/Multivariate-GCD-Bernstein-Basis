function opt_col = GetOptimalColumn(fxy_matrix_n,gxy_matrix_n,t1,t2)
%% Find Optimal column for removal from S_{t_{1},t_{2}}
% Given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized

% Inputs

% fxy_matrix : Matrix of coefficients of polynomial f(x,y)
%
% gxy_matrix : Matrix of coefficients of polynomial g(x,y)
%
% t1         : Degree of GCD d(x,y) with respect to x
%
% t2         : Degree of GCD d(x,y) with respect to y


% Get the degree of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix_n);

St1t2 = BuildDTQ(fxy_matrix_n,gxy_matrix_n,t1,t2);

%% Find Optimal column for removal from St
% given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized

% Get the number
[~,ncolsSt1t2] = size(St1t2);

% QR Decomposition of the Sylvester Matrix S_{k}
[Qk,Rk] = qr(St1t2);

n = n1+n2;

% For each column of S_{t_{1},t_{2}}, move the column to the RHS, to obtain
% Ax=b, where A consists of the remaining columns of S_{t_{1},t_{2}} and b
% is the removed column c_{k}

for k = 1 : 1 : ncolsSt1t2
    
    % Get column for removal
    ck = St1t2(:,k);
    
    % Perform QR delete to remove k column from QR decomposition of 
    % S_{t_{1},t_{2}}
    [Q,~] = qrdelete(Qk,Rk,k);
    
    cd = Q'*ck;
    
    d = cd(n+1:end,:);
    
    % Get Residuals
    residuals_QR(k) = norm(d);
    
end

%Obtain the column for which the residual is minimal.
[~,opt_col] = min(log10(residuals_QR));

% Print out optimal column for removal.
fprintf([mfilename ' : ' sprintf('Optimal column for removal is : %i \n',opt_col)]);

end