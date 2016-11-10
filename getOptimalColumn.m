function idx_col = GetOptimalColumn(Sk)
%% Find Optimal column for removal from S_{t_{1},t_{2}}
% Given that t1 and t2 have been calculated build the sylvester matrix and
% find the optimal column such that a residual is minimized
%
% % Inputs
%
% Sk : Sylvester subresultant matrix S_{k}(f,g)
%
% % Outputs
%
% idx_col : Index of optimal column


% Get the number of columns in S_{k}(f,g)
[~,nCols_Sk] = size(Sk);

% Initialise vector 
vResiduals = zeros(nCols_Sk,1);

for k = 1 : 1 : nCols_Sk
    
    % Get column for removal
    ck = Sk(:,k);
    Ak = Sk;
    Ak(:,k) = [];
    
    % Get the vector x.
    x = SolveAx_b(Ak,ck);
    
    % Get Residuals
    vResiduals(k) = norm(ck - Ak*x);
    
end

%Obtain the column for which the residual is minimal.
[~,idx_col] = min(log10(vResiduals));

% Print out optimal column for removal.
fprintf([mfilename ' : ' sprintf('Optimal column for removal is : %i \n',idx_col)]);

end