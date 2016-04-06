function [t1,t2,lambda,mu,alpha,th1,th2] = GetGCDDegree_Relative(fxy,gxy,...
    m,n,t,lambda,mu,alpha,th1,th2)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two 
% polynomials f(x,y) and g(x,y)
%
%   Inputs.
%
%   fxy : Coefficient matrix of polynomial f(x,y)
%
%   gxy : Coefficient matrix of polynomial g(x,y)
%
%   m : Total degree of polynomial f(x,y)
%
%   n : Total degree of polynomial g(x,y)
%
%   t : Total degree of GCD d(x,y)
%
%   lambda : Geometric mean of entries in first partition of Syvlester
%   matrix
%
%   mu : Geometric mean of Entries in second partition of Sylvester matrix
%
%   alpha : Optimal value of alpha
%
%   th1 : Optimal value of theta_{1}
%
%   th2 : Optimal value of theta_{2}


global PLOT_GRAPHS

% Get the degree structure of polynomial f(x,y)
[m1,m2] = GetDegree(fxy);

% Get the degree structure of polynomial g(x,y)
[n1,n2] = GetDegree(gxy);

% Produce the set of all possible t1 and t2 values
method = 'All';

switch method
    case 'All'  % Use all possible (t1,t2) combinations
        % The total of t1+t2 must be between t and 2t
        mat = [];
        for t1 = min(m1,n1):-1:0;
            for t2 = min(m2,n2):-1:0;
                
                new_row = [t1 t2];
                mat = [mat ; new_row];
                
            end
        end
        
    case 'Refined' % Use only a subset of possible (t1,t2) combinations.
        % The total of t1+t2 must be between t and 2t
        mat = [];
        
        for t1 = t:-1:0;
            for t2 = t:-1:0;
                
                condA = n1-t1 + n2 -t2 >= n-t;
                condB = m1-t1 + m2 -t2 >= m-t;
                condC = t1 <= n1 && t1 <= m1;
                condD = t2 <= n2 && t2 <= m2;
                condE = n1 - t1 + n2 - t2 <= 2*(n-t);
                condF = m1 - t1 + m2 - t2 <= 2*(m-t);
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
[nRowsMyMat,~] = size(mat);

% for each row in the constructed matrix of (k1,k2) pairs
data = [];


% for every row in the matrix
for i = 1:1:nRowsMyMat
    
    k1 = mat(i,1);
    k2 = mat(i,2);
    
    % Preprocessing
    [lambda, mu, alpha, th1,th2] = Preprocess(fxy,gxy,k1,k2);
    
    % Divide f(x) by geometric mean
    fxy_matrix_n = fxy ./lambda;
    
    % Divide g(x) by geometric mean
    gxy_matrix_n = gxy ./mu;
    
    % Get f(w,w) from f(x,y)
    fww_matrix = GetWithThetas(fxy_matrix_n,th1,th2);
    
    % Get g(w,w) from g(x,y)
    gww_matrix = GetWithThetas(gxy_matrix_n,th1,th2);
    
    % Build the k1,k2 subresultant matrix
    Sk1k2 = ...
        BuildDTQ(fww_matrix,alpha.*gww_matrix,k1,k2);
    
    % Get the minimum singular value
    min_sing_val = min(svd(Sk1k2));
    
    
    try
        data = [data ; k1 k2 log10(min_sing_val) alpha th1 th2 lambda mu ];
        z(k1+1,k2+1) = log10(min_sing_val);
    catch
        data = [data ; k1 k2 0];
        z(k1+1,k2+1) = 0;
    end
    
    
end



max_k1 = max(data(:,1));
max_k2 = max(data(:,2));


switch method
    case 'All'
        y = 0:min(m2,n2)+1;
        x = 0:min(m1,n1)+1;
        
    case 'Refined'
        x = 0:max_k1;
        y = 0:max_k2;
        
end

[nRowsMyMat,nColsMyMat] = size(z);
z2 = zeros(nRowsMyMat+1,nColsMyMat+1);
z2(1:nRowsMyMat,1:nColsMyMat) = z;
z = z2;

%%
% Plot 3d surface
switch PLOT_GRAPHS
    case 'y'
        [x,y] = meshgrid(x,y);
        figure_name = sprintf('%s - Surface',mfilename);
        figure('name',figure_name)
        hold on
        s1 = surf(x,y,z');
        xlabel('t_{1}')
        ylabel('t_{2}')
        xlim([0,max_k1+3])
        ylim([0,max_k2+3])
        
        xlabel('t_{1}')
        ylabel('t_{2}')
        hold off
    case 'n'
    otherwise
        error('Error : plot_graphs is either (y) or (n)')
end

%%
switch PLOT_GRAPHS
    case 'y'
        % Plot 3d data points
        figure_name = sprintf('%s - Minimum Singular Values', mfilename);
        figure('name',figure_name)
        hold on
        title('Minimum Singular Values in S_{t_{1},t_{2}}')
        xlabel('t_{1}')
        ylabel('t_{2}')
        zlabel('log_{10} Min Sing Value')
        scatter3(data(:,1),data(:,2),data(:,3),'filled')
        grid('on')
        hold off
    case 'n'
    otherwise
        error('Error: plot_graphs is either y or n')
end
%%
[nRowsMyMat,nColsMyMat] = size(z);
% take from (0,0)
zz = z(1:1:nRowsMyMat,1:1:nColsMyMat);

% get the second row to the end
delta_x = [zeros(1,nColsMyMat) ; z(2:1:nRowsMyMat,1:1:nColsMyMat)];

% get the second col to the end
delta_y = [zeros(nRowsMyMat,1) z(1:1:nRowsMyMat,2:1:nColsMyMat)];

% Get change in x component (k1) and a row of zeros
delta_z_x = [diff(z,1) ; zeros(1,nColsMyMat)];

% Get change in y component (k2) and a zero col

delta_z_y = [diff(z,1,2) zeros(nRowsMyMat,1)];

delta_zz = delta_z_x + delta_z_y;

criterion = max(log10(abs(delta_zz(:))));

global THRESHOLD
if criterion < THRESHOLD
    fprintf('Value below threshold \n')
    
    t1 = min(m1,n1);
    t2 = min(m2,n2);
    
    display(t1)
    display(t2)
    
    return
end

% Get the [i,j] entry which has maximum change to [i+1,j] and [i,j+1]
[num idx] = max(delta_zz(:));
[x y] = ind2sub(size(delta_zz),idx);
% Set the values t1 and t2
t1 = x-1;
t2 = y-1;


% % If only one row is returned, then take t1 and t2
[nRowsMyMat,~] = size(data);
if nRowsMyMat == 1
    t1 = data(1,1);
    t2 = data(1,2);
    return;
end

% fprintf('The set of pairs t_{1} and t_{2} are given by: \n')
% mymat

% Get the position of the maximum change in values of min singular value
% for each t1 + t2 = tot

switch PLOT_GRAPHS
    case 'y'
        figure_name = sprintf('%s - data',mfilename);
        figure('name',figure_name)
        plot(log10(data(:,2)));
        hold off
    case 'n'
    otherwise
        error('plot_graphs is either y or n')
end


fprintf('\n')
fprintf('Degree t1 : %i \n',t1);
fprintf('Degree t2 : %i \n',t2);
fprintf('\n')
fprintf('----------------------------------------------------------\n')


end




