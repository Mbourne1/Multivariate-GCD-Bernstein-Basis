function [] = o_roots_mymethod(fxy_matrix)
% o_roots_mymethod(fxy_matrix)
%
% Get the roots of the input polynomial f(x,y)
%
% Inputs.
%
%
% fxy_matrix : Matrix of coefficients of polynomial f(x,y)


%
%               Root Finding Algorithm

% Obtain the series q_{x}{i} by series of GCD calculations

% Set the first entry of q to be the input polynomial f(x,y)
qx{1} = fxy_matrix;

% Get dimensions of polynomial f(x,y)
[r,c] = size(fxy_matrix);
m1 = r-1;
m2 = c-1;

% Set the degree structure of qx{1}
deg1_qx{1} = m1;
deg2_qx{1} = m2;

% Set the iteration number
ite_num = 1;

% Set the iteration condition to true
iteration_condition = true;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while iteration_condition
    str1 = iptnum2ordinal(ite_num);
    str2 = sprintf('GCD calculation with respect to x. \n',ite_num);
    fprintf([str1 str2])
            
    % Get the degree structure of f(x,y)
    [r,c] = size(qx{ite_num});
    m1 = r - 1;
    m2 = c - 1;
    m = m1 + m2;
    
    % Differentiate f(x,y) with respect to x to obtain g(x,y)
    qx_der = Bern_Differentiate_wrt_x(qx{ite_num});
    n1 = m1-1;
    n2 = m2;
    n = n1+n2;
    
    fprintf('Input Polynomials for GCD calculation : %i \n',ite_num)
    
    % If the degree of g(x,y) with respect to x is 0 , then g(x,y) is a
    % scalar and set GCD to be g(x,y)
    if n1 == 0
        % GCD is only a scalar with respect to x so set equal to g(x,y).
        [~,~,dxy_matrix,t1,t2] = o1(qx{ite_num},qx_der,m,n);
    else
        % Get the GCD of f(x,y) and g(x,y)
        [~,~,dxy_matrix,t1,t2] = o1(qx{ite_num},qx_der,m,n); 
    end
    
    % increment the iteration number
    ite_num = ite_num + 1;
    
    % Assign the GCD to be the newest member of the q array.
    qx{ite_num} = dxy_matrix;
    deg1_qx{ite_num} = t1;
    deg2_qx{ite_num} = t2;
    
    % set m to be equal to t, ready for the next iteration
    m1 = t1;
    m2 = t2;
    
    % If the degree of the GCD with respect to x is zero, then terminate 
    % iterations
    if t1 == 0
        iteration_condition = false;
    else
        iteration_condition = true;
    end
end


%% Obtain the series h_{x}{i} by series of deconvolutions on q_{x}{i}

% Get number of elements in the series of polynomials q_{i}
[~,num_entries_qx] = size(qx);

% Pre assign the CellArray h_{x} to have one less element than the
% CellArray q_{x}
hx = cell(1,num_entries_qx-1);
deg1_hx = cell(1,num_entries_qx-1);
deg2_hx = cell(1,num_entries_qx-1);

% For each pair of consecutive polynomials in qx perform deconvolution
for i = 1:1:num_entries_qx-1
    
    % Get the series h_{x,i}
    hx{i} = Bernstein_Deconvolve_Bivariate(qx{i},qx{i+1});
    
    % Set the degree of h with respect to x
    deg1_hx{i} = deg1_qx{i} - deg1_qx{i+1};
    % Set the degree of h with respect to y
    deg2_hx{i} = deg2_qx{i} - deg2_qx{i+1};
    
end


% Get the number of entries in the array of h_{x}.
[~,num_entries_hx] = size(hx);

%% Obtain the series w_{x}{i} by series of deconvolutions on h_{x}{i}

% Pre assign the CellArray w_{x}
num_entries_wx = num_entries_hx -1;

wx = cell(1,num_entries_wx);
deg1_wx = cell(1,num_entries_wx);
deg2_wx = cell(1,num_entries_wx);

% For each pair, perform a deconvolution to obtain w_{x}
if num_entries_hx > 1
    for i = 1:1: num_entries_hx - 1
         
        wx{i} = Bernstein_Deconvolve_Bivariate(hx{i},hx{i+1});
        % Set the degree of w with respect to x
        deg1_wx{i} = deg1_hx{i} - deg1_hx{i+1};
        % Set the degree of w with respect to y
        deg2_wx{i} = deg2_hx{i} - deg2_hx{i+1};
    end
    
    % Include the final element of hx in wx
    wx{i+1} = hx{i+1};
    % Set its degree with respect to x
    deg1_wx{i+1} = deg1_hx{i+1};
    % Set its degree with respect to y
    deg2_wx{i+1} = deg2_hx{i+1};
    
else
    
    % number of elements in the array h_{i} is equal to one
    wx{1} = hx{1};
    % Set the degree with respect to x
    deg1_wx{1} = deg1_hx{1};
    % Set the degree with respect to y
    deg2_wx{1} = deg2_hx{1};
end

fprintf('################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obtain series of q_{y}{i} by series of GCD calculations

% Set the first entry of q to be the input polynomial f(x,y)
qy{1} = fxy_matrix;
deg1_qy{1} = m1;
deg2_qy{1} = m2;

% Get dimensions of polynomial q_{y}(x,y)
[r,c] = size(qy{1});
m1 = r-1;
m2 = c-1;

% Set the iteration number
ite_num = 1;

% If f(x,y) is just a scalar with respect to y then dont perform GCD
% calculation
if m2 < 1
    iteration_condition = false;
else
    iteration_condition = true;
end

%

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while iteration_condition
    
    str1 = iptnum2ordinal(ite_num);
    str2 = sprintf(' GCD calculation with respect to y. %i \n',ite_num);
    fprintf([str1 str2])
    
    % set dxy to be the input of the next loop
    fxy_matrix = qy{ite_num};
    
    [r,c] = size(qy{ite_num});
    m1 = r - 1;
    m2 = c - 1;
    m = m1 + m2;
    
    % Differentiate f(x,y) with respect to x to obtain g(x,y)
    qy_der = Bern_Differentiate_wrt_y(fxy_matrix);
    [r,c] = size(qy_der);
    n1 = r - 1;
    n2 = c - 1;
    n = n1 + n2;
       
    % if g(x,y) is a scalar with respect to y   
    if n2 == 0
        % GCD is only a scalar with respect to y so set equal to g(x,y).
        [uxy,vxy,dxy_matrix,t1,t2] = o1(qy{ite_num},qy_der,m,n);
    else
        % Get the GCD of f(x,y) and g(x,y)
        [uxy,vxy,dxy_matrix,t1,t2] = o1(qy{ite_num},qy_der,m,n);
    end
    
    
    
    % increment the iteration number
    ite_num = ite_num + 1;
    
    % Assign the GCD to be the newest member of the q array.
    qy{ite_num} = dxy_matrix;
    deg1_qy{ite_num} = t1;
    deg2_qy{ite_num} = t2;

    
    % Check the iteration condition.
    if t2 == 0
        iteration_condition = false;
    else
        iteration_condition = true;
    end
end

%% Obtain the series h_{y}{i}

% get number of elements in the series of polynomials q_{i}
[~,num_entries_qy] = size(qy);

num_entries_hy = num_entries_qy - 1;

% Preassign CellArray for hy
hy = cell(1,num_entries_hy)
deg1_hy = cell(1,num_entries_hy);
deg2_hy = cell(1,num_entries_hy);

% For every q_{y,i}, deconvolve with q_{y,i+1} to obtain h_{y,i}
for i = 1:1:num_entries_qy-1
    
    % Perform Deconvolution to obtain h_{y,i}
    hy{i} = Bernstein_Deconvolve_Bivariate(qy{i},qy{i+1});
    
    % Get the degree of h_{y,i} with respect to x
    deg1_hy{i} = deg1_qy{i} - deg1_qy{i+1};
    % Get the degree of h_{y,i} with respect to y
    deg2_hy{i} = deg2_qy{i} - deg2_qy{i+1};
   
    
end



%% obtain the series w_{i} for the

if exist('hy') ~= 0
    [~,c] = size(hy);
    
    if c > 1
        for i = 1:1:c-1
            wy{i}    = Bernstein_Deconvolve_Bivariate(hy{i},hy{i+1});
            deg1_wy{i} = deg1_hy{i} - deg1_hy{i+1};
            deg2_wy{i} = deg2_hy{i} - deg2_hy{i+1};
        end
        
        wy{i+1} = hy{i+1};
        deg1_wy{i+1} = deg1_hy{i+1};
        deg2_wy{i+1} = deg2_hy{i+1};
    else
        wy{1} = hy{1};
        deg1_wy{1} = deg1_hy{1};
        deg2_wy{1} = deg2_hy{1};
    end
end


% Print out the set of w_{i} in term of x
if exist('wx') ~=0
    fprintf('----------------------------------------------------------------\n')
    fprintf('Printing the roots with respect to x\n')
    if exist('wx') ~=0
        % For every entry in w_{x}
        for i = 1:1:num_entries_wx
            str1 = iptnum2ordinal(i);
            str2 = sprintf(' entry in w_{x}. \n');
            fprintf(['The ' str1 str2])
            wxi = cell2mat(wx(i));
            wxi./wxi(1,1)
        end
    end
end

if exist('wy') ~= 0
    fprintf('----------------------------------------------------------------\n')
    fprintf('Printing the roots with respect to y\n')
    [r,c] = size(wy);
    for i = 1:1:c
        str1 = iptnum2ordinal(i);
        str2 = sprintf(' entry in w_{y}. \n');
        fprintf(['The ' str1 str2])
        wyi = cell2mat(wy(i));
        wyi./wyi(1,1)
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Perform a series of GCD calculations on the w_{x,i}s

% get number of w_{x,i}
for i = 1:1:num_entries_wx
    
    % Get the dimensions of the polynomial the ith polynomial in w_{x}
    [rows,cols] = size(wx{i});
    m1 = rows - 1;
    m2 = cols - 1;
    
    % If the polynomial has a y component (is Bivariate), then must 
    % deconvolve with the equivalent polynomial w_{y}
    
    if m2 > 0
        
        [uxy_calc_matrix, vxy_calc_matrix, dxy_calc_matrix,t1,t2]  =...
            o1(wx{i},wy{i},deg1_wx{i},deg1_wy{i});
        
        % Overwrite wx and wy with new values
        % Assign the GCD to the non-Separable part
        wxy{i} = dxy_calc_matrix;
        % Divide w_{x} by the GCD to obtain w_{x} without y component
        wx{i} = Bernstein_Deconvolve_Bivariate(wx{i},dxy_calc_matrix);
        % Divide w_{y} by the GCD to obtain w_{y} without x component
        wy{i} = Bernstein_Deconvolve_Bivariate(wy{i},dxy_calc_matrix);
        
    end
    
end

%%
% For each w_{x,i} get the root
% get the polynomial, whose roots have multiplicty i, in bernstein
% form, where coefficients are in terms of (1-y)^{m-i}y^{i}.
fprintf('----------------------------------------------------------------\n')
fprintf('Roots and Multiplicities of f(x,y) with respect to x \n')
[r,c] = size(wx);
for i=1:1:c
    
    % Get the factor
    factor = wx{i};
    [r2,c2] = size(factor);
    
    if (r2 == 2 )
        aw = wx{i};
        
        % Normalise the factors polynomial coefficients
        aw = aw./aw(1);
        
        % Convert to power form, so that coefficients are in terms of y^{i}
        % rather than (1-y)^{m-i}y^{i}.
        a_pwr = [aw(1,:) ; aw(2,:)-aw(1,:)];
        
        % Obtain the root in terms of y, and set multiplicity to one.
        fprintf('The root of multiplicity %i is given by:\n', i)
        [-a_pwr(1,:)./a_pwr(2,:) a_pwr(2,:)./a_pwr(2,:)]
    end
    
end
fprintf('----------------------------------------------------------------\n')
%%
% For each w_{y,i} get the root
% get the polynomial, whose roots have multiplicty i, in bernstein
% form, where coefficients are in terms of (1-y)^{m-i}y^{i}.
if exist('wy') ~= 0
    fprintf('Roots and Multiplicities of f(x,y) with respect to y \n')
    [r,c] = size(wy);
    for i=1:1:c

        % Get the factor
        factor = wy{i};
        [r2,c2] = size(factor);

        if (c2 == 2 )
            aw = wy{i};

            % Normalise the factors polynomial coefficients
            aw = aw./aw(1);

            % Convert to power form, so that coefficients are in terms of y^{i}
            % rather than (1-y)^{m-i}y^{i}.
            a_pwr = [aw(:,1) ; aw(:,2)-aw(:,1)];

            % Obtain the root in terms of y, and set multiplicity to one.
            fprintf('The root of multiplicity %i is given by:\n', i)
            [-a_pwr(1,:)./a_pwr(2,:) a_pwr(2,:)./a_pwr(2,:)]
        end

    end
end
fprintf('----------------------------------------------------------------\n')

%%
% For each w_{y,i} get the root
% get the polynomial, whose roots have multiplicty i, in bernstein
% form, where coefficients are in terms of (1-y)^{m-i}y^{i}.
if exist('wxy') ~=0
    fprintf('Bivariate Factors of f(x,y) \n')
    [r,c] = size(wxy);
    for i=1:1:c

        fprintf('Factor of multiplicity %i is given by:\n', i)
        % Get the factor
        wxy{i}



    end
end

end