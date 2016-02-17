function P = BuildP(m1,m2,n1,n2,theta1,theta2,alpha,t1,t2,opt_col)

num_cols_T1 = (n1 - t1 + 1) * (n2 - t2 + 1);

% Get the number of coefficients in the polynomial f(x,y)
num_coeff_f = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
num_coeff_g = (n1+1) * (n2+1);

if opt_col <= num_cols_T1
    
    % Optimal column in first partition
    fprintf('Optimal column in First partition\n')
    
    % % Build the Matrix P
    % Build the matrix P1
    P1 = BuildP1(m1,m2,n1,n2,theta1,theta2,opt_col,t1,t2);
    
    % Build the matrix P2
    rows = (m1+n1-t1+1)*(m2+n2-t2+1);
    P2 = zeros(rows,num_coeff_g);
 

else
    % Optimal column in second partition
    fprintf('Optimal column in second partition\n')
    % % Build the Matrix P
    
    % Build the matrix P1
    rows = (m1+n1-t1+1)*(m2+n2-t2+1);
    P1 = zeros(rows,num_coeff_f);
    
    % Build the matrix P2
    % Get the position of the optimal column with respect to T(g)
    opt_col_rel = opt_col - num_cols_T1;
    P2 = BuildP1(n1,n2,m1,m2,theta1,theta2,opt_col_rel,t1,t2);
    

end

P =  [P1 alpha.*P2];

end