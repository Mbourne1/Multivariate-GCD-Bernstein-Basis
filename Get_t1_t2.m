function [t1,t2,lambda,mu,opt_alpha, opt_theta_1,opt_theta_2] = Get_t1_t2(fxy_matrix,gxy_matrix,...
    m,n,t,...
    lambda,mu,...
    opt_alpha,opt_theta_1,opt_theta_2)

global bool_preproc

[r,c] = size(fxy_matrix);
m1 = r-1;
m2 = c-1;

[r,c] = size(gxy_matrix);
n1 = r-1;
n2 = c-1;


%%
% given that we know t
% it is possible that t1 and t2 take any values such that
%   1.) t1 + t2 >= t
%   2.) t1 + t2 <= 2t
%   3.) t1 <= t
%   4.) t2 <= t
%   2.) t1 <= m1
%   3.) t1 <= n1
%   4.) t2 <= m2
%   5.) t2 <= n2

%    t1 <= min(m1,n1)
%    t2 <= min(m2,n2)


% Produce the set of all possible t1 and t2 values
method = 'All';

switch method
    case 'All'
        % The total of t1+t2 must be between t and 2t
        mat = [];
        for t1 = min(m1,n1):-1:0;
            for t2 = min(m2,n2):-1:0;
                
                new_row = [t1 t2];
                mat = [mat ; new_row];
                
            end
        end
        
    case 'Refined'
        % The total of t1+t2 must be between t and 2t
        mat = [];
        
        for t1 = t:-1:0;
            for t2 = t:-1:0;
                
                condA = n1-t1 + n2 -t2 >= n-t;
                condB = m1-t1 + m2 -t2 >= m-t;
                condC = t1 <= n1 && t1 <= m1;
                condD = t2 <= n2 && t2 <= m2;
                condE = n1 - t1 + n2 - t2 <= 2*(n-t);
                condF = m1 - t1 + n1 - t1 <= 2*(m-t);
                condG = n1 - t1 <= n-t;
                condH = n2 - t2 <= n-t;
                condI = m1 - t1 <= m-t;
                condJ = m2 - t2 <= m-t;
                
                
                if condA && condB && condC && condD && condE && condF...
                        && condG && condH && condI && condJ
                    new_row = [t1 t2];
                    mat = [mat ; new_row];
                end
                
            end
        end
end

% Remove duplicate rows
mat = unique(mat,'rows');
[r,~] = size(mat);

% for each row in the constructed matrix of (k1,k2) pairs
mymat = [];


% for every row in the matrix
for i = 1:1:r
    
    k1 = mat(i,1);
    k2 = mat(i,2);
    
    %% Apply preprocessing
    switch bool_preproc
        case 'y'
            
            % Preproecessor One - Normalise by geometric mean
            [lambda, mu] = getGeometricMean(fxy_matrix,gxy_matrix,k1,k2);
            
            % Normalise f(x,y)
            fxy_matrix_n = fxy_matrix./lambda;
            
            % Normalise g(x,y)
            gxy_matrix_n = gxy_matrix./mu;
            
            % Preprocessor Two and Three - LinProg to obtain optimal values
            % of alpha, theta_1 and theta_2
            
            % Get the maximum and minimum entries of f(x,y) in the
            % Sylvester matrix S(f,g)
            [max_mtrx_f, min_mtrx_f] = GetMaxMin(fxy_matrix_n,n1,n2,k1,k2);
            
            % Get the maximum and minimum entries of g(x,y) in the
            % Sylvester matrix S(f,g)
            [max_mtrx_g, min_mtrx_g] = GetMaxMin(gxy_matrix_n,m1,m2,k1,k2);
            
            % Get optimal values of alpha and theta
            [opt_alpha, opt_theta_1, opt_theta_2] = OptimalAlphaTheta(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g);
            
        case 'n'
            fxy_matrix_n = fxy_matrix;
            gxy_matrix_n = gxy_matrix;
            opt_alpha = 1;
            opt_theta_1 = 1;
            opt_theta_2 = 1;
            lambda = 1;
            mu = 1;
            
    end
    
    % Build the k1,k2 subresultant
    Sk1k2 = ...
        BuildSubresultant(fxy_matrix_n,gxy_matrix_n,k1,k2,opt_alpha,opt_theta_1,opt_theta_2);
    
    % Get the minimum singular value
    min_sing_val = min(svd(Sk1k2));
    
    
    try
        mymat = [mymat ; k1 k2 log10(min_sing_val) opt_alpha opt_theta_1 opt_theta_2 lambda mu ];
        z(k1+1,k2+1) = log10(min_sing_val);
    catch
        mymat = [mymat ; k1 k2 0];
        z(k1+1,k2+1) = 0;
    end
    
    
end



max_k1 = max(mymat(:,1));
max_k2 = max(mymat(:,2));


switch method
    case 'All'
        y = 0:min(m2,n2)+1;
        x = 0:min(m1,n1)+1;
        
    case 'Refined'
        x = 0:max_k1;
        y = 0:max_k2;
        
end

[r,c] = size(z);
z2 = zeros(r+1,c+1);
z2(1:r,1:c) = z;
z = z2;

%%
% Plot 3d surface
[x,y] = meshgrid(x,y);
figure('name','surface')
hold on
s1 = surf(x,y,z');
xlabel('t_{1}')
ylabel('t_{2}')
xlim([0,max_k1+3])
ylim([0,max_k2+3])
alpha(s1,0.5)
xlabel('t_{1}')
ylabel('t_{2}')
hold off


%%
% Plot 3d data points
figure('name','3dplot')
hold on
title('Minimum Singular Values in S_{t_{1},t_{2}}')
xlabel('t_{1}')
ylabel('t_{2}')


scatter3(mymat(:,1),mymat(:,2),mymat(:,3),'filled')
grid('on')
hold off

%%
[r,c] = size(z);
% take from (0,0)
zz = z(1:1:r,1:1:c);

% get the second row to the end
delta_x = [zeros(1,c) ; z(2:1:r,1:1:c)];

% get the second col to the end
delta_y = [zeros(r,1) z(1:1:r,2:1:c)];

% Get change in x component (k1) and a row of zeros
delta_z_x = [diff(z,1) ; zeros(1,c)];

% Get change in y component (k2) and a zero col

delta_z_y = [diff(z,1,2) zeros(r,1)];

delta_zz = delta_z_x + delta_z_y;

criterion = max(log10(abs(delta_zz(:))));

global threshold
if criterion < threshold
    fprintf('Value below threshold \n')
    t1 = min(m1,n1)
    t2 = min(m2,n2)
    
    return
end

% Get the [i,j] entry which has maximum change to [i+1,j] and [i,j+1]
[num idx] = max(delta_zz(:));
[x y] = ind2sub(size(delta_zz),idx);
% Set the values t1 and t2
t1 = x-1;
t2 = y-1;


% % If only one row is returned, then take t1 and t2
[r,~] = size(mymat);
if r == 1
    t1 = mymat(1,1);
    t2 = mymat(1,2);
    return;
end

% fprintf('The set of pairs t_{1} and t_{2} are given by: \n')
% mymat

% Get the position of the maximum change in values of min singular value
% for each t1 + t2 = tot


figure('name','Separate Degree Calculation')
plot(log10(mymat(:,2)));
hold off








