function [t1,t2,opt_alpha, opt_theta_1_matrix,opt_theta_2_matrix] = GetDegree(fxy_matrix,gxy_matrix,m,n)

% Get the degree of the two bivariate bernstein basis polynomials.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs.


% fxy_matrix :

% gxy_matrix :

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Outputs.


% t1 :  Degree of GCD with respect to x

% t2 :  Degree of GCD with respect to y

% opt_theta_1 : Optimal value of theta with respect to x =
%               \theta_{1}\omega{1}

% opt_theta_2 : Optimal value of theta with respect to y =
%               \theta_{2}\omega_{2}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global bool_preproc


% Get degrees of fxy 
[rows,cols] = size(fxy_matrix);
% Get degree of f in terms of x
m1 = rows - 1;
% Get degree of f in terms of y
m2 = cols - 1;


% Get degree of gxy
[rows,cols] = size(gxy_matrix);
% Get degree of g in terms of x
n1 = rows - 1;
% Get degree of g in terms of y
n2 = cols - 1;

% Initialise some empty vectors for storing during the loop
sing_val_vec = zeros(min(m,n),1);
t_val_vec1 = zeros(min(m,n),2);
t_val_vec2 = zeros(min(m,n),2);
t_val_vec3 = zeros(min(m,n),2);

% initialise some useful vectors
Data_RowNorm    = [];
Data_DiagNorm   = [];

% let k represent the 
for k = 1:1:min(m,n)
    min_sing_val_vec = [];
    temp_vec = [];
    count = 1;
    ratio_maxmin_diag_vec = [];
    ratio_maxmin_rowsum_vec = [];
    
    
    for k1 = 0:1:k
        k2 = k-k1;
        if (k1 > m1 || k1 > n1) || (k2 > m2 || k2 > n2)
        else
            
            % STAGE 1 - Preprocess by geometric mean
            
            switch bool_preproc
                case 'y'
                    [GM_f,GM_g] = getGeometricMean(fxy_matrix,gxy_matrix,k1,k2);
                case 'n'
                    GM_f = 1;
                    GM_g = 1;
                    
            end
            
            % Obtain normalised fxy and gxy matrices where
            % binomials are included.
            %fxy_matrix_bi_n = fxy_matrix_bi./GM_f;
            %gxy_matrix_bi_n = gxy_matrix_bi./GM_g;
            
            % Obtain normalised fxy and gxy matrices.
            fxy_matrix_n = fxy_matrix./GM_f;
            gxy_matrix_n = gxy_matrix./GM_g;
            
            
            switch bool_preproc
                case 'y'
                    % Get the maximum and minimum entries of each
                    % coefficient of f, in the Sylvester matrix
                    [max_mtrx_f,min_mtrx_f] = GetMaxMin(fxy_matrix_n,n1,n2,k1,k2);
                    
                    % Get the maximum and minimum entries of each
                    % coefficient of g, in the Sylvester matrix
                    [max_mtrx_g,min_mtrx_g] = GetMaxMin(gxy_matrix_n,n1,n2,k1,k2);
                    
                    % Get an optimal value of alpha alone
                    %opt_alpha = OptimalAlpha(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g)
                    
                    % Get optimal values of theta1 and theta2 alone
                    %[opt_theta_1, opt_theta_2] = OptimalTheta(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g);
                    
                    % Get optimal values of theta1 and theta2 and alpha
                    [opt_alpha, opt_theta_1, opt_theta_2] = OptimalAlphaTheta(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g)
                    
                    % Add the optimal values of theta 1 and theta 2 to a
                    % matrix defined by k1,k2
                    
                case 'n' % Exclude preprocessors
                    % set values of theta to 1.
                    opt_theta_1 = 1;
                    opt_theta_2 = 1;
                    opt_alpha = 1;
            end
            
            opt_theta_1_matrix(k1+1,k2+1) = opt_theta_1;
            opt_theta_2_matrix(k1+1,k2+1) = opt_theta_2;
            opt_alpha_matrix(k1+1,k2+1) = opt_alpha;
            
            % Build the subresultant s_{k_{1},k_{2}}
            Sk = BuildSubresultant(...
                fxy_matrix_n, gxy_matrix_n,...
                k1, k2,...
                opt_alpha,...
                opt_theta_1,...
                opt_theta_2);
            
            
            % Edit 24/07/2015
            % Using QR Decomposition of the sylvester matrix
            [~,R] = qr(Sk);
            
            % Take absolute values.
            R = abs(R);
            
            % Get number of rows in R1
            [R1_rows,~] = size(diag(R));
            
            % Obtain R1 the top square of the R matrix.
            R1 = R(1:R1_rows,1:R1_rows);
            
            % Get Norms of each row in the matrix R1
            R1_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);
            
            % Get ONLY the diagonal elements and normalise them.
            R1_DiagNorm = diag(R1)./norm(diag(R1));
            
            % Scatter Plot Data
            ks = count.*ones(size(R1_RowNorm));
            ns = 1:1:size(R1_RowNorm,1);
            
            % Form a triple of [ks, the value of QR_RowNorm, and the index of the value of
            % the row of R1 corresponding to QR_RowNorm].
            % EG.
            %  [1   0.015  1
            %   1   0.156  2
            %   2 ...]
            X = [ks R1_RowNorm ns'];
            Data_RowNorm = [Data_RowNorm; X];
            
            X2 = [ks R1_DiagNorm ns'];
            Data_DiagNorm = [Data_DiagNorm;X2];
            
            % Get ratio of max diag elem of S_{k} to min diag elemente of S_{k}
            ratio_maxmin_diag_vec = [ratio_maxmin_diag_vec ; max(diag(R1))./min(diag(R1))];
            ratio_maxmin_rowsum_vec = [ratio_maxmin_rowsum_vec ; max(R1_RowNorm)./min(R1_RowNorm)];
            
            % End Edit 24/07/2015
            %             figure()
            %             hold on
            %             plot(ratio_maxmin_diag_vec)
            %             hold off
            
            % add the minimum singular value S_{k1,k2} to the vector of
            % singular values for all k = k1+k2
            min_sing_val_vec(count) = min(svd(Sk));
            temp_vec(count,:) = [k1 k2];
            count = count + 1;
        end
    end
    try
    % of all the diagonal entries in R1 for the given k1+k2 = k, get the
    % largest ratio
    [ratio_maxmin_diag(k),index1] = max(ratio_maxmin_diag_vec);
    t_val_vec1(k,:) = temp_vec(index1,:);
    
    [ratio_maxmin_rowsum(k),index2] = max(ratio_maxmin_rowsum_vec);
    t_val_vec2(k,:) = temp_vec(index2,:);
    
    % of all the minimal singular values for the given k1+k2 = k, stored in
    % the vector min_sing_val_vec, get the minimal singular value and its
    % index.
    [sing_val_vec(k),index] = min(min_sing_val_vec);
    
    % get the values of k1 and k2 which gave the minimal value
    t_val_vec3(k,:) = temp_vec(index,:);
    catch 
    end
end

opt_theta_1_matrix;
opt_theta_2_matrix;

% plot all the largest ratios for k = 1,...,min(m,n)
figure()
hold on
plot(log10(ratio_maxmin_diag));
hold off

figure()
hold on
plot(log10(ratio_maxmin_rowsum));
hold off

% plot all the minimum singular values for k = 1,...,min(m,n)
figure()
hold on
plot(log10(sing_val_vec))
hold off


diff(sing_val_vec);

[~,maxindex] = max(diff(log10(sing_val_vec)));
degree_calc = maxindex;

[tval] = t_val_vec3(degree_calc,:);
t1 = tval(1);
t2 = tval(2);


fprintf('The Calculated Degree of the GCD is given by \n %i \n \n',degree_calc)
fprintf('t1 = %i\n',t1)
fprintf('t2 = %i\n',t2)

% Get the optimal values of alpha and theta used to build the subresultant
% S_{t_{1},t_{2}}
opt_theta1 = opt_theta_1_matrix(t1+1,t2+1)
opt_theta2 = opt_theta_2_matrix(t1+1,t2+1)

end