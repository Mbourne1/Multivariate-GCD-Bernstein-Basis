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

% %
%               Root Finding Algorithm

% Obtain the series q_{x}{i} by series of GCD calculations

% Set the first entry of q to be the input polynomial f(x,y)
fx{1} = fxy_matrix;

% Get dimensions of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get total degree of f(x,y)
m = m1+m2;

% Set the degree structure of qx{1}
vDeg1_fx(1) = m1;
vDeg2_fx(1) = m2;
vDegt_fx(1) = m;

% Set the iteration number
ite = 1;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDeg1_fx(ite) > 0
    
    fprintf('GCD Calculation Loop iteration = %i \n', ite );
    fprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite);   
    
    % Differentiate f(x,y) with respect to x to obtain g(x,y)
    gxy = Differentiate_wrt_x(fx{ite});
    
    m = vDegt_fx(ite);
    n  = m-1;
    
    fprintf('Input Polynomials for GCD calculation : %i \n',ite)
    
    
    if ite > 1
        lower_lim = vDegt_fx(ite)-d(ite-1);
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    fprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim);
    fprintf('Maximum degree of f_{%i}: %i \n', ite+1, upper_lim);
    
    % GCD is only a scalar with respect to x so set equal to g(x,y).
    [fx{ite},~,fx{ite+1},t,t1,t2] = o_gcd_mymethod(fx{ite},gxy,m,n,[lower_lim,upper_lim]);

    % Set the degree of q{i} with respect to x
    vDeg1_fx(ite+1) = t1;
    
    % Set the degree of q{i} with respect to y
    vDeg2_fx(ite+1) = t2;
    
    % Set the total degree of q{i} 
    vDegt_fx(ite+1) = t;
    
    % Get number of distinct roots of f(ite)
    d(ite) = vDegt_fx(ite) - vDegt_fx(ite+1);
    
    fprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_fx(ite+1))
    fprintf('Number of distinct roots in f_{%i} : %i \n',ite,d(ite))
    fprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_fx(ite+1))
    
    % increment the iteration number
    ite = ite + 1;
end


%% Obtain the series h_{x}{i} by series of deconvolutions on q_{x}{i}

% Get number of elements in the series of polynomials q_{i}
[~,num_entries_qx] = size(fx);

% Pre assign the CellArray h_{x} to have one less element than the
% CellArray q_{x}
hx = cell(1,num_entries_qx-1);

% Initialise vectors to store the degrees of h(x).
vDeg1_hx = zeros(1,num_entries_qx-1);
vDeg2_hx = zeros(1,num_entries_qx-1);
vDegt_hx = zeros(1,num_entries_qx-1);

% For each pair of consecutive polynomials in qx perform deconvolution
for i = 1:1:num_entries_qx-1
    
    % Get the series h_{x,i}
    hx{i} = Deconvolve_Bivariate(fx{i},fx{i+1});
    
    % Set the degree of h with respect to x
    vDeg1_hx(i) = vDeg1_fx(i) - vDeg1_fx(i+1);
    
    % Set the degree of h with respect to y
    vDeg2_hx(i) = vDeg2_fx(i) - vDeg2_fx(i+1);
    
    % Set the total degree of h 
    vDegt_hx(i) = vDegt_fx(i) - vDegt_fx(i+1);
    
end


% Get the number of entries in the array of h_{x}.
[~,num_entries_hx] = size(hx);

%% Obtain the series w_{x}{i} by series of deconvolutions on h_{x}{i}

% Pre assign the CellArray w_{x}
num_entries_wx = num_entries_hx -1;

wx = cell(1,num_entries_wx);

% Initialise vectors to store the degrees of w(x).
vDeg1_wx = zeros(1,num_entries_wx);
vDeg2_wx = zeros(1,num_entries_wx);
vDegt_wx = zeros(1,num_entries_wx);

% For each pair, perform a deconvolution to obtain w_{x}
if num_entries_hx > 1
    for i = 1:1: num_entries_hx - 1
         
        wx{i} = Deconvolve_Bivariate(hx{i},hx{i+1});
        
        % Set the degree of w with respect to x
        vDeg1_wx(i) = vDeg1_hx(i) - vDeg1_hx(i+1);
        
        % Set the degree of w with respect to y
        vDeg2_wx(i) = vDeg2_hx(i) - vDeg2_hx(i+1);
        
        % Set the total degree
        vDegt_wx(i) = vDegt_hx(i) - vDegt_hx(i+1);
    end
    
    % Include the final element of hx in wx
    wx{i+1} = hx{i+1};
    
    % Set its degree with respect to x
    vDeg1_wx(i+1) = vDeg1_hx(i+1);
    
    % Set its degree with respect to y
    vDeg2_wx(i+1) = vDeg2_hx(i+1);
    
    % Set the total degree
    vDegt_wx(i+1) = vDegt_hx(i+1);
    
else
    
    % number of elements in the array h_{i} is equal to one
    wx{1} = hx{1};
    
    % Set the degree with respect to x
    vDeg1_wx(1) = vDeg1_hx(1);
    
    % Set the degree with respect to y
    vDeg2_wx(1) = vDeg2_hx(1);
    
    % Set the total degree
    vDegt_wx(1) = vDegt_hx(1);
end

fprintf('################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obtain series of q_{y}{i} by series of GCD calculations

% Set the first entry of q to be the input polynomial f(x,y)
fy{1} = fxy_matrix;

[m1,m2] = GetDegree(fxy_matrix);
m = m1+m2;

% Get degree of polynomial f(x,y) with respect to x
vDeg1_fy(1) = m1;

% Get degree of polynomial f(x,y) with respect to y
vDeg2_fy(1) = m2;

% Get total degree of polynomial f(x,y)
vDegt_fy(1) = m;

% Set the iteration number
ite = 1;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.

while vDegt_fy(ite) > 0
    
   
    % set dxy to be the input of the next loop
    fxy_matrix = fy{ite};
    
    % Get the degree structure
    [m1,m2] = GetDegree(fy{ite});
    
    m = m1 + m2;
    
    % Differentiate f(x,y) with respect to x to obtain g(x,y)
    qy_der = Differentiate_wrt_y(fxy_matrix);
    
    % Get degree of g(x,y)
    n = m-1;
    
    if ite > 1
        lower_lim = vDegt_fy(ite)-d(ite-1);
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    fprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim);
    fprintf('Maximum degree of f_{%i}: %i \n', ite+1, upper_lim);
    
    % GCD is only a scalar with respect to y so set equal to g(x,y).
    [fy{ite},~,fy{ite+1},t,t1,t2] = o_gcd_mymethod(fy{ite},qy_der,m,n,[lower_lim,upper_lim]);
      
    %
    vDeg1_fy(ite+1) = t1;
    
    %
    vDeg2_fy(ite+1) = t2;
    
    %
    vDegt_fy(ite+1) = t;
    
    % increment the iteration number
    ite = ite + 1;
    
end

%% Obtain the series h_{y}{i}

% get number of elements in the series of polynomials q_{i}
[~,num_entries_qy] = size(fy);

num_entries_hy = num_entries_qy - 1;

% Preassign CellArray for hy
hy = cell(1,num_entries_hy);

%
vDeg1_hy = zeros(1,num_entries_hy);

%
vDeg2_hy = zeros(1,num_entries_hy);

% For every q_{y,i}, deconvolve with q_{y,i+1} to obtain h_{y,i}
for i = 1:1:num_entries_qy-1
    
    % Perform Deconvolution to obtain h_{y,i}
    hy{i} = Deconvolve_Bivariate(fy{i},fy{i+1});
    
    % Get the degree of h_{y,i} with respect to x
    vDeg1_hy(i) = vDeg1_fy(i) - vDeg1_fy(i+1);
    
    % Get the degree of h_{y,i} with respect to y
    vDeg2_hy(i) = vDeg2_fy(i) - vDeg2_fy(i+1);
   
    % Get the total degree of h_{y,i}
    vDegt_hy(1) = vDegt_fy(i) - vDegt_fy(i+1);
end



%% obtain the series w_{i}(y)

if exist('hy') ~= 0
    [~,c] = size(hy);
    
    if c > 1
        for i = 1:1:c-1
            wy{i}    = Deconvolve_Bivariate(hy{i},hy{i+1});
            vDeg1_wy(i) = vDeg1_hy(i) - vDeg1_hy(i+1);
            vDeg2_wy(i) = vDeg2_hy(i) - vDeg2_hy(i+1);
            vDegt_wy(i) = vDegt_hy(i) - vDegt_hy(i+1);
        end
        
        wy{i+1} = hy{i+1};
        vDeg1_wy(i+1) = vDeg1_hy(i+1);
        vDeg2_wy(i+1) = vDeg2_hy(i+1);
        vDegt_wy(i+1) = vDegt_hy(i+1);
    else
        wy{1} = hy{1};
        vDeg1_wy(1) = vDeg1_hy(1);
        vDeg2_wy(1) = vDeg2_hy(1);
        vDegt_wy(1) = vDegt_hy(1);
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
    [~,nCols_wy] = size(wy);
    for i = 1:1:nCols_wy
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
    [~,m2] = GetDegree(wx{i});
        
    % If the polynomial has a y component (is Bivariate), then must 
    % deconvolve with the equivalent polynomial w_{y}
    
    if m2 > 0
        
        [uxy_calc_matrix, vxy_calc_matrix, dxy_calc_matrix,t,t1,t2]  =...
            o1(wx{i},wy{i},vDegt_wx(i),vDegt_wy(i));
        
        % Overwrite wx and wy with new values
        % Assign the GCD to the non-Separable part
        wxy{i} = dxy_calc_matrix;
        
        % Divide w_{x} by the GCD to obtain w_{x} without y component
        wx{i} = Deconvolve_Bivariate(wx{i},dxy_calc_matrix);
        
        % Divide w_{y} by the GCD to obtain w_{y} without x component
        wy{i} = Deconvolve_Bivariate(wy{i},dxy_calc_matrix);
        
    end
    
end

%%
% For each w_{x,i} get the root
% get the polynomial, whose roots have multiplicty i, in bernstein
% form, where coefficients are in terms of (1-y)^{m-i}y^{i}.
fprintf('----------------------------------------------------------------\n')
fprintf('Roots and Multiplicities of f(x,y) with respect to x \n')
[~,c] = size(wx);
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