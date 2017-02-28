function [wy, vDegt_wy] = o_roots_mymethod_y(fxy, M)
% Set the first entry of q to be the input polynomial f(x,y)
fy{1} = fxy;

% Get degree of polynomial f(x,y)
[m1,m2] = GetDegree_Bivariate(fxy);

% Get degree of polynomial f(x,y) with respect to x
vDeg1_fy(1) = m1;

% Get degree of polynomial f(x,y) with respect to y
vDeg2_fy(1) = m2;

% Get total degree of polynomial f(x,y)
vDegt_fy(1) = M;

% Set the iteration number
ite = 1;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.

while vDeg2_fy(ite) > 0
    
    if (vDeg2_fy(ite) == 1)
        % The derivative is a constant
        
        % The GCD is a constant
        fy{ite+1} = Differentiate_wrt_y(fy{ite});
        uy{ite+1} = Deconvolve_Bivariate(fy{ite},fy{ite+1});
        
        % Get degree of d(x,y) with respect to x
        vDeg1_fy(ite+1) = 0;
        
        % Get degree of d(x,y) with respect to y
        vDeg2_fy(ite+1) = vDeg2_fy(ite);
        
        % Get total degree of d(x,y)
        vDegt_fy(ite+1) = 0;
        break;
    end
        
    fprintf([mfilename ' : ' sprintf('GCD Calculation Loop iteration = %i \n', ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite)]);
    
    % Differentiate f(x,y) with respect to x to obtain g(x,y)
    gxy = Differentiate_wrt_y(fy{ite});
    
    % Get total degree of f(x,y)
    m = vDegt_fy(ite);
    
    % Get total degree of g(x,y)
    n = m-1;
    
    % Set upper and lower bounds of total degree.
    if ite > 1
        lower_lim = vDegt_fy(ite)-d(ite-1);
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n', ite+1, upper_lim)]);
    
    % GCD is only a scalar with respect to y so set equal to g(x,y).
    [fy{ite},gxy,fy{ite+1},uxy,vxy,t,t1,t2] = o_gcd_mymethod_2Polys(fy{ite}, gxy, m, n, [lower_lim,upper_lim]);
      
    % Set the degree of q{i} with respect to x
    vDeg1_fy(ite+1) = t1;
    
    % Set the degree of q{i} with respect to y
    vDeg2_fy(ite+1) = t2;
    
    % Set the total degree of q{i}
    vDegt_fy(ite+1) = t;
    
    % Get number of distinct roots of f(ite)
    d(ite) = vDegt_fy(ite) - vDegt_fy(ite+1);
    
    fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_fy(ite+1))]);
    fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,d(ite))]);
    fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_fy(ite+1))]);
    LineBreakLarge();
    
    % increment the iteration number
    ite = ite + 1;
    
end

%% Obtain the series h_{y}{i}

% Get number of elements in the series of polynomials q_{i}
[~,nEntries_qy] = size(fy);

nEntries_hy = nEntries_qy - 1;

% Preassign CellArray for h_{y} to have one less element that the 
% CellArray q_{y}
hy = cell(1,nEntries_hy);

% Initialise vectors to store the degrees of h(y)
vDeg1_hy = zeros(1,nEntries_hy);
vDeg2_hy = zeros(1,nEntries_hy);
vDegt_hy = zeros(1,nEntries_hy);

% For every q_{y,i}, deconvolve with q_{y,i+1} to obtain h_{y,i}
for i = 1:1:nEntries_hy
    
    % Perform Deconvolution to obtain h_{y,i}
    hy{i} = Deconvolve_Bivariate(fy{i},fy{i+1});
    
    % Get the degree of h_{y,i} with respect to x
    vDeg1_hy(i) = vDeg1_fy(i) - vDeg1_fy(i+1);
    
    % Get the degree of h_{y,i} with respect to y
    vDeg2_hy(i) = vDeg2_fy(i) - vDeg2_fy(i+1);
   
    % Get the total degree of h_{y,i}
    vDegt_hy(1) = vDegt_fy(i) - vDegt_fy(i+1);
end

% Get the number of entires in the array h_{y}
[~,nEntries_hy] = size(hy);

%% obtain the series w_{i}(y) by series of deconvolutions on h_{y}{i}

nEntries_wy = nEntries_hy -1 ;

wy = cell(1,nEntries_wy);

vDeg1_wy = zeros(1,nEntries_wy);
vDeg2_wy = zeros(1,nEntries_wy);
vDegt_wy = zeros(1,nEntries_wy);

% for each pair, perform a deconvolution to obtain w_{x}
if nEntries_hy > 1
    for i = 1:1:nEntries_hy - 1

        wy{i} = Deconvolve_Bivariate(hy{i},hy{i+1});
        
        % Set the degree of w with respect to x
        vDeg1_wy(i) = vDeg1_hy(i) - vDeg1_hy(i+1);
        
        % Set the degree of w with respect to y
        vDeg2_wy(i) = vDeg2_hy(i) - vDeg2_hy(i+1);
        
        % Set the total degree
        vDegt_wy(i) = vDegt_hy(i) - vDegt_hy(i+1);
    end
         % Include the final element of hx in wx
        wy{i+1} = hy{i+1};

        % Set its degree with respect to x
        vDeg1_wy(i+1) = vDeg1_hy(i+1);

        % Set its degree with respect to y
        vDeg2_wy(i+1) = vDeg2_hy(i+1);

        % Set the total degree
        vDegt_wy(i+1) = vDegt_hy(i+1);
        
else
    
    % number of elements in the array h_{i} is equal to one
    wy{1} = hy{1};
    
    % Set the degree with respect to x
    vDeg1_wy(1) = vDeg1_hy(1);
    
    % Set the degree with respect to y
    vDeg2_wy(1) = vDeg2_hy(1);
    
    % Set the total degree
    vDegt_wy(1) = vDegt_hy(1);
end

end