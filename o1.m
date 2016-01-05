function [uxy_calc_matrix, vxy_calc_matrix, dxy_calc_matrix,t1,t2] = o1(fxy_matrix,gxy_matrix,...
    m,n)
% Given two input polynomials, calculate the GCD and its degree structure
%% Get the degree of the GCD

% Get degree calculation method
degree_calc_method = '1';
switch degree_calc_method
    % Get total degree
    case '1'
        % Get Degree by first finding the total degree, then obtain t1 and t2
        
        % Get total degreee
        [t, opt_theta_1, opt_theta_2] = Get_t(fxy_matrix, gxy_matrix,m,n);
        
        fprintf('---------------------------------------------------------\n')
        fprintf('\n')
        fprintf('Total Degree of GCD as calculated by GetDegree_Total() : \n')
        fprintf('Tota Degree = %i',t)
        fprintf('\n')
        fprintf('---------------------------------------------------------\n')
        
        
        % Get degree t1 and t2
        [t1,t2,lambda,mu,opt_alpha, opt_theta_1,opt_theta_2] = Get_t1_t2(fxy_matrix,gxy_matrix,m,n,t);
        
        fprintf('\n')
        fprintf('Degree t1 : %i \n',t1);
        fprintf('Degree t2 : %i \n',t2);
        fprintf('\n')
        fprintf('----------------------------------------------------------\n')
    otherwise
        error('no other degree calc method is defined')
end

%% Get Optimal column for removal from S_{t_{1},t_{2}}

opt_col = getOptimalColumn(fxy_matrix,gxy_matrix,t1,t2,lambda,mu,opt_alpha,opt_theta_1,opt_theta_2);


%% Perform iterative inprovements in SNTLN
global bool_SNTLN

switch bool_SNTLN
    case 'y'
        % Apply SNTLN improvements
        [ fxy_output,gxy_output,alpha_output,theta1_output,theta2_output,X_output] = ...
            SNTLN( fxy_matrix,gxy_matrix,...
            opt_alpha, opt_theta_1, opt_theta_2,...
            t1,t2,...
            lambda,mu,...
            opt_col);
        
        fprintf('Input Polynomial f(x,y)')
        
        
    case 'n'
        % Dont Apply SNTLN improvements
        
    otherwise
        error('bool_SNTLN is either y or n')
end

%% Get Quotients u(x,y) and v(x,y)

[uxy_calc_matrix, vxy_calc_matrix, lambda,mu, opt_alpha, opt_theta_1, opt_theta_2] ...
    = GetQuotients(fxy_matrix, gxy_matrix,t1,t2);





%% Get the GCD
% % Get d(x,y) from the polynomials u(x,y) and v(x,y).
dxy_calc_matrix = GetGCD_Coefficients(uxy_calc_matrix,vxy_calc_matrix,...
    fxy_matrix,gxy_matrix,...
    t1,t2,...
    lambda,mu,...
    opt_alpha,opt_theta_1,opt_theta_2);



end




