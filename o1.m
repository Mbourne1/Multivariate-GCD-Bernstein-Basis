function [uxy_calc_matrix, vxy_calc_matrix, dxy_calc_mtrx] = o1(fxy_matrix_working,gxy_matrix_working,...
                                    m,n,...
                                    t_exct,t1_exct,t2_exct,...
                                    uxy_mtrx_exct,vxy_mtrx_exct,dxy_mtrx_exct,...
                                    fxy_matrix_exact,gxy_matrix_exact)                                
%% Get Initial data
[r,c] = size(fxy_matrix_exact);
m1 = r - 1;
m2 = c - 1;

[r,c] = size(gxy_matrix_exact);
n1 = r - 1;
n2 = c - 1;

                                
                                
%% Get the degree of the GCD

% Get degree calculation method
degree_calc_method = input('Degree Calculation method \n(1) By total degree \n(2) by t1 t2 \n Method : ','s');

switch degree_calc_method
% Get total degree
    case '1'
        %% Get Degree by first finding the total degree, then obtain t1 and t2
        
        % Get total degreee
        [t,opt_theta_1, opt_theta_2] = GetDegree_Total(fxy_matrix_working, gxy_matrix_working,m,n);
        
        fprintf('Total Degree of GCD as calculated by GetDegree_Total() : \n')
        fprintf('tdeg = %i',t)
        fprintf('\n')
        fprintf('Optimal values of theta: \n')
        fprintf('Theta 1 : %0.5e \n',opt_theta_1)
        fprintf('Theta 2 : %0.5e \n',opt_theta_2)
        fprintf('\n')
        % Get degree t1 and t2
        GetDegree_Total_Separate(fxy_matrix_working,gxy_matrix_working,m,n,t)
        
        % Return Optimal alpha and theta
        
        
    case '2'
%%
        % Get the degrees t1 and t2
        % Also return the corresponding values of \theta_{1} and \theta_{2}
        [t1,t2, opt_theta_1_mtrx, opt_theta_2_mtrx] = GetDegree(fxy_matrix_working,gxy_matrix_working,...
            m,n);
        % Get the optimal values of theta for the subresualtant S_{t1,t2}
        opt_theta1 = opt_theta_1_mtrx(t1+1,t2+1);
        opt_theta2 = opt_theta_2_mtrx(t1+1,t2+1);
end
%%
% TEMPORARY 
t1 = t1_exct;
t2 = t2_exct;
t = t_exct;
%%


%% Get Quotients

[uxy_calc_matrix, vxy_calc_matrix, opt_theta_1, opt_theta_2] = GetQuotients(fxy_matrix_working, gxy_matrix_working,t1,t2);

fprintf('----------------------------------------------------------------')
fprintf('\n')
fprintf('Calculated and exact values of u(x,y):')
uxy_calc_matrix ./ uxy_calc_matrix(1,1)
uxy_mtrx_exct./uxy_mtrx_exct(1,1)
fprintf('----------------------------------------------------------------')
fprintf('\n')
fprintf('Calculated and exact values of v(x,y): ')

vxy_calc_matrix ./ vxy_calc_matrix(1,1)
vxy_mtrx_exct./vxy_mtrx_exct(1,1)
fprintf('----------------------------------------------------------------')
fprintf('\n')


%% Get the GCD
% % Get d from the polynomials u and v.


dxy_calc_mtrx = GetGCD_Coefficients(uxy_calc_matrix,vxy_calc_matrix,...
    fxy_matrix_working,gxy_matrix_working,...
    t1,t2,...
    opt_theta_1,opt_theta_2);

fprintf('----------------------------------------------------------------')
fprintf('\n')
fprintf('Calculated and Exact values of d(x,y):')
d1 = dxy_calc_mtrx./dxy_calc_mtrx(1,1)
d2 = dxy_mtrx_exct./dxy_mtrx_exct(1,1)

%% Get the GCD as a vector

% Build d calc as a vector
dxy_vec_calc = getAsVector(dxy_calc_mtrx)

% Build d exact as a vector
dxy_vec_exact = getAsVector(dxy_mtrx_exct)


fprintf('Compare Exact coefficients with computed coefficients')
[(dxy_vec_calc./dxy_vec_calc(1)) ./ (dxy_vec_exact./dxy_vec_exact(1))  ]


end




