function [wx,vDegt_wx] = o_roots_mymethod_x(fxy_matrix,M)
% o_roots_mymethod_x(fxy_matrix,M)
%
% Inputs
%
% fxy_matrix : Coefficients of polynomial f(x,y)
%
% M : Degree of f(x,y)


% Set the first entry of q to be the input polynomial f(x,y)
fx{1} = fxy_matrix;

% Get degree of polynomial f(x,y)
%[m1,m2] = GetDegree(fxy_matrix);

% Get degree of polynomial f(x,y) with respect to x
%vDeg1_fx(1) = m1;

% Get degree of polynomial f(x,y) with respect to y
%vDeg2_fx(1) = m2;

% Get total degree of polynomial f(x,y)
vDegt_fx(1) = M;

% Set the iteration number
ite = 1;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDeg1_fx(ite) > 0
    
    if (vDeg1_fx(ite) == 1)
        % The derivative is a constant
        
        % The GCD is a constant
        fx{ite+1} = Differentiate_wrt_x(fx{ite});
        ux{ite+1} = Deconvolve_Bivariate(fx{ite},fx{ite+1});
        
        % Get degree of d(x,y) with respect to x
        %vDeg1_fx(ite+1) = 0;
        
        % Get degree of d(x,y) with respect to y
        %vDeg2_fx(ite+1) = vDeg2_fx(ite);
        
        % Get total degree of d(x,y)
        vDegt_fx(ite+1) = 0;
        break;
    end
    
    fprintf([mfilename ' : ' sprintf('GCD Calculation Loop iteration = %i \n', ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite)]);
    
    % Differentiate f(x,y) with respect to x to obtain g(x,y)
    gxy = Differentiate_wrt_x(fx{ite});
    
    % Get total degree of f(x,y)
    m = vDegt_fx(ite);
    
    % Get total degree of g(x,y)
    n  = m-1;
   
    % Set upper and lower bounds of total degree.
    if ite > 1
        lower_lim = vDegt_fx(ite)-d(ite-1);
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n', ite+1, upper_lim)]);
    
    % GCD is only a scalar with respect to x so set equal to g(x,y).
    [fx{ite},gxy,fx{ite+1},uxy,vxy,t,t1,t2] = o_gcd_mymethod(fx{ite},gxy,m,n,[lower_lim,upper_lim]);
    
    % Set the degree of q{i} with respect to x
    %vDeg1_fx(ite+1) = t1;
    
    % Set the degree of q{i} with respect to y
    %vDeg2_fx(ite+1) = t2;
    
    % Set the total degree of q{i}
    vDegt_fx(ite+1) = t;
    
    % Get number of distinct roots of f(ite)
    d(ite) = vDegt_fx(ite) - vDegt_fx(ite+1);
    
    fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_fx(ite+1))]);
    fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,d(ite))]);
    fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_fx(ite+1))]);
    LineBreakLarge();
    
    % Increment the iteration number
    ite = ite + 1;
end


%% Obtain the series h_{x}{i} by series of deconvolutions on q_{x}{i}

% Get number of elements in the series of polynomials q_{i}
[~,nEntries_qx] = size(fx);

nEntries_hx = nEntries_qx - 1;

% Pre assign the CellArray h_{x} to have one less element than the
% CellArray q_{x}
hx = cell(1,nEntries_hx);

% Initialise vectors to store the degrees of h(x).
%vDeg1_hx = zeros(1,nEntries_hx);
%vDeg2_hx = zeros(1,nEntries_hx);
vDegt_hx = zeros(1,nEntries_hx);

% For each pair of consecutive polynomials in qx perform deconvolution
for i = 1:1:nEntries_hx
    
    % Get the series h_{x,i}
    hx{i} = Deconvolve_Bivariate(fx{i},fx{i+1});
    
    % Set the degree of h with respect to x
    %vDeg1_hx(i) = vDeg1_fx(i) - vDeg1_fx(i+1);
    
    % Set the degree of h with respect to y
    %vDeg2_hx(i) = vDeg2_fx(i) - vDeg2_fx(i+1);
    
    % Set the total degree of h
    vDegt_hx(i) = vDegt_fx(i) - vDegt_fx(i+1);
    
end


% Get the number of entries in the array of h_{x}.
[~,nEntries_hx] = size(hx);

%% Obtain the series w_{x}{i} by series of deconvolutions on h_{x}{i}

% Pre assign the CellArray w_{x}
nEntries_wx = nEntries_hx -1;

wx = cell(1,nEntries_wx);

% Initialise vectors to store the degrees of w(x).
%vDeg1_wx = zeros(1,nEntries_wx);
%vDeg2_wx = zeros(1,nEntries_wx);
vDegt_wx = zeros(1,nEntries_wx);

% For each pair, perform a deconvolution to obtain w_{x}
if nEntries_hx > 1
    for i = 1:1: nEntries_hx - 1
        
        wx{i} = Deconvolve_Bivariate(hx{i},hx{i+1});
        
        % Set the degree of w with respect to x
        %vDeg1_wx(i) = vDeg1_hx(i) - vDeg1_hx(i+1);
        
        % Set the degree of w with respect to y
        %vDeg2_wx(i) = vDeg2_hx(i) - vDeg2_hx(i+1);
        
        % Set the total degree
        vDegt_wx(i) = vDegt_hx(i) - vDegt_hx(i+1);
    end
    
    % Include the final element of hx in wx
    wx{i+1} = hx{i+1};
    
    % Set its degree with respect to x
    %vDeg1_wx(i+1) = vDeg1_hx(i+1);
    
    % Set its degree with respect to y
    %vDeg2_wx(i+1) = vDeg2_hx(i+1);
    
    % Set the total degree
    vDegt_wx(i+1) = vDegt_hx(i+1);
    
else
    
    % number of elements in the array h_{i} is equal to one
    wx{1} = hx{1};
    
    % Set the degree with respect to x
    %vDeg1_wx(1) = vDeg1_hx(1);
    
    % Set the degree with respect to y
    %vDeg2_wx(1) = vDeg2_hx(1);
    
    % Set the total degree
    vDegt_wx(1) = vDegt_hx(1);
end

end
