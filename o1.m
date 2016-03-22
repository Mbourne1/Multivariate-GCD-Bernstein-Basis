function [uxy_calc_matrix, vxy_calc_matrix, dxy_calc_matrix,t1,t2] = o1(fxy_matrix,gxy_matrix,...
    m,n)
% o1(fxy_matrix,gxy_matrix,m,n)
%
% Given two input polynomials, calculate the GCD and its degree structure
%
% Inputs.
%
%
% fxy_matrix : Matrix of coefficients of f(x,y)
%
% gxy_matrix : Matrix of coefficients of g(x,y)
%
% m : Total degree of polynomial f(x,y).
%
% n : Total degree of polynomial g(x,y).



% Get the degree of the GCD

degree_calc_method = 'respective';

% Get Degree by first finding the total degree, then obtain t1 and t2

% Get total degreee
[t, opt_theta_1_tot, opt_theta_2_tot] = Get_t(fxy_matrix, gxy_matrix,m,n);

% Get degree t1 and t2
[t1,t2,lambda,mu,alpha, th1,th2] = Get_t1_t2(fxy_matrix,gxy_matrix,m,n,t);

% Get Optimal column for removal from S_{t_{1},t_{2}}
opt_col = GetOptimalColumn(fxy_matrix,gxy_matrix,t1,t2,lambda,mu,alpha,th1,th2);


%% Perform iterative inprovements in SNTLN
global LOW_RANK_APPROXIMATION_METHOD

switch LOW_RANK_APPROXIMATION_METHOD
    case 'Standard SNTLN'
        % Apply SNTLN improvements
        [ fxy_matrix,gxy_matrix,alpha_output,theta1_output,theta2_output,X_output] = ...
            SNTLN( fxy_matrix,gxy_matrix, alpha, th1, th2,t1,t2, lambda,mu, opt_col);
        
        fprintf('Input Polynomial f(x,y)')
        
    case 'Standard STLN'
        %% Preprocessing

        % Normalise polynomial f(x,y) by geometric means
        fxy_matrix_n = fxy_matrix./lambda;
        gxy_matrix_n = gxy_matrix./mu;

        % Obtain polynomials in Modified Bernstein Basis, using initial values of
        % alpha and theta.

        % Multiply the rows of fxy_matrix by theta1, and multiply the cols of
        % fxy_matrix by theta2
        fww_matrix = GetWithThetas(fxy_matrix_n,th1,th2);

        % Multiply the rows of gxy_matrix by theta1, and multiply the cols of
        % gxy_matrix by theta2
        gww_matrix = GetWithThetas(gxy_matrix_n,th1,th2);
        [fww_matrix, gww_matrix, X_output] = STLN(fww_matrix, alpha.*gww_matrix,t1,t2,opt_col)
        
        fxy_matrix = GetWithoutThetas(fww_matrix,th1,th2);
        gxy_matrix = GetWithoutThetas(gww_matrix,th1,th2) ./ alpha;
        
        display('')
        
    case 'None'
        % Dont Apply SNTLN improvements
        
    otherwise
        error('bool_SNTLN is either (Standard SNTLN) or (None)')
end

%% Get Quotients u(x,y) and v(x,y)
% Calc method is either total or respective

switch degree_calc_method
    case 'respective'
        [uxy_calc_matrix, vxy_calc_matrix, lambda,mu, alpha, th1, th2] ...
            = GetQuotients(fxy_matrix, gxy_matrix,t1,t2);
    case 'total'
        [uxy_calc_matrix, vxy_calc_matrix, lambda,mu, alpha, th1, th2] ...
            = GetQuotients_total(fxy_matrix, gxy_matrix,m,n,t);
    otherwise
        
end




%% Get the GCD
% % Get d(x,y) from the polynomials u(x,y) and v(x,y).

switch degree_calc_method
    case 'respective'
        dxy_calc_matrix = GetGCD_Coefficients(uxy_calc_matrix,vxy_calc_matrix,...
            fxy_matrix,gxy_matrix,...
            t1,t2,...
            lambda,mu,...
            alpha,th1,th2);
    case 'total'
        dxy_calc_matrix = GetGCD_Coefficients_total(uxy_calc_matrix,vxy_calc_matrix,...
            fxy_matrix,gxy_matrix,...
            m,n,t,...
            lambda,mu,...
            alpha,th1,th2);
    otherwise
        error('error')
end


end




