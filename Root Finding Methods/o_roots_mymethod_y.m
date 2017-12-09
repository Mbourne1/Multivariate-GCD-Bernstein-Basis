function [arr_wxy] = o_roots_mymethod_y(fxy)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% % Outputs
%
% arr_wxy : (Array of Matrices) Contain matrices of coefficients of
% polynomials w_{i}(x,y)
%
% vDegree_x_wxy : (Vector) Degree of polynomials w_{i}(x,y) with respect to
% x
%
% vDegree_y_wxy : (Vector) Degree of polynomials w_{i}(x,y) with respect to
% y
%
% % Outputs
%
% arr_wxy : (Array of Matrices)
%
% vDegree_x_wxy :
%
% vDegree_y_wxy :


% Get array of polynomials f_{i}(x,y)
arr_fxy = GetArray_fxy(fxy);

% Get array of polynomials h_{i}(x,y)
arr_hxy = GetArray_hxy(hxy);

% Get array of polynomials w_{i}(x,y)
arr_wxy = GetArray_wxy(wxy);









end


function arr_fxy = GetArray_fxy(fxy)
%
% % Inputs
% 
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% % Outputs
%
% arr_fxy : (Array of Matrices) 


% Set the first entry of q to be the input polynomial f(x,y)
arr_fxy{1} = fxy;

% Get degree of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degree of polynomial f(x,y) with respect to x
vDegree_x_fxy(1) = m1;

% Get degree of polynomial f(x,y) with respect to y
vDegree_y_fxy(1) = m2;

% Set the iteration number
ite = 1;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.

while vDegree_y_fxy(ite) > 0
    
    if (vDegree_y_fxy(ite) == 1)
        % The derivative is a constant
        
        % The GCD is a constant
        arr_fxy{ite + 1} = Differentiate_wrt_y(arr_fxy{ite});
        arr_uxy{ite + 1} = Deconvolve_Bivariate(arr_fxy{ite},arr_fxy{ite+1});
        
        % Get degree of d(x,y) with respect to x
        vDegree_x_fxy(ite + 1) = 0;
        
        % Get degree of d(x,y) with respect to y
        vDegree_y_fxy(ite + 1) = vDegree_y_fxy(ite);
        
        break;
    end
    
    fprintf([mfilename ' : ' sprintf('GCD Calculation Loop iteration = %i \n', ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite)]);
    
    % Differentiate f(x,y) with respect to x to obtain g(x,y)
    gxy = Differentiate_wrt_y(arr_fxy{ite});
    
    % Get degree of f(x,y) and g(x,y)
    [m1, m2] = GetDegree_Bivariate(fxy);
    [n1, n2] = GetDegree_Bivariate(gxy);
    
    % Set upper and lower bounds of total degree.
    if ite > 1
        
        lowerLimit_t1 = vDegree_x_fxy(ite) - dx(ite-1);
        upperLimit_t1 = min(m1, n1);
        
        lowerLimit_t2 = vDegree_y_fxy(ite) - dy(ite-1);
        upperLimit_t2 = min(m2, n2);
        
        
        
    else
        lowerLimit_t1 = 1;
        upperLimit_t1 = min(m1, n1);
        
        lowerLimit_t2 = 1;
        upperLimit_t2 = min(m2, n2);
        
        
    end
    
    limits_t1 = [lowerLimit_t1, upperLimit_t1];
    limits_t2 = [lowerLimit_t2, upperLimit_t2];
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i} with respect to x: %i \n', ite+1, lowerLimit_t1)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i} with respect to y: %i \n', ite+1, upperLimit_t1)]);
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i} with respect to y: %i \n', ite+1, lowerLimit_t2)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i} with respect to y: %i \n', ite+1, upperLimit_t2)]);
    
    % GCD is only a scalar with respect to y so set equal to g(x,y).
    [arr_fxy{ite}, gxy, arr_fxy{ite+1}, uxy, vxy, t1, t2] = o_gcd_mymethod_Bivariate_2Polys(arr_fxy{ite}, gxy, limits_t1, limits_t2);
    
    % Set the degree of q{i} with respect to x
    vDegree_x_fxy(ite+1) = t1;
    
    % Set the degree of q{i} with respect to y
    vDegree_y_fxy(ite+1) = t2;
    
    dx(ite) = vDegree_x_fxy(ite) - vDegree_x_fxy(ite+1);
    dy(ite) = vDegree_y_fxy(ite) - vDegree_y_fxy(ite+1);
    
    %fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegree_t_fxy(ite+1))]);
    fprintf([mfilename ' : ' sprintf('The degree of (GCD(f_{%i},f_{%i}) with respect to x is : %i \n',ite, ite, vDegree_x_fxy(ite+1))]);
    fprintf([mfilename ' : ' sprintf('The degree of (GCD(f_{%i},f_{%i}) with respect to y is : %i \n',ite, ite, vDegree_y_fxy(ite+1))]);
    %fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,d(ite))]);
    %fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite + 1, vDegree_t_fxy(ite+1))]);
    LineBreakLarge();
    
    % increment the iteration number
    ite = ite + 1;
    
end


end





function arr_hxy = GetArray_hxy(arr_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices)
%
% % Outputs
%
% arr_hxy : (Array of Matrices)


% % Obtain the series h_{y}{i}

% Get number of elements in the series of polynomials q_{i}
nPolys_fxy = length(arr_fxy);

nPolys_hxy = nPolys_fxy - 1;

% Preassign CellArray for h_{y} to have one less element that the
% CellArray q_{y}
arr_hxy = cell(1, nPolys_hxy);

% Initialise vectors to store the degrees of h(y)
%vDegree_x_hxy = zeros(1, nPolys_hxy);
%vDegree_y_hxy = zeros(1, nPolys_hxy);

% For every q_{y,i}, deconvolve with q_{y,i+1} to obtain h_{y,i}
for i = 1:1:nPolys_hxy
    
    % Perform Deconvolution to obtain h_{y,i}
    arr_hxy{i} = Deconvolve_Bivariate(arr_fxy{i}, arr_fxy{i+1});
    
    % Get the degree of h_{y,i} with respect to x
    %vDegree_x_hxy(i) = vDegree_x_fxy(i) - vDegree_x_fxy(i+1);
    
    % Get the degree of h_{y,i} with respect to y
    %vDegree_y_hxy(i) = vDegree_y_fxy(i) - vDegree_y_fxy(i+1);
    
end


end



function arr_wxy = GetArray_wxy(arr_hxy)
%
% % Inputs
%
% arr_hxy : (Array of Matrices)
%
% % Outputs
%
% arr_wxy : (Array of Matrices)




% % obtain the series w_{i}(y) by series of deconvolutions on h_{y}{i}

nEntries_wxy = nEntries_hxy -1 ;

arr_wxy = cell(1,nEntries_wxy);

%vDegree_x_wxy = zeros(1,nEntries_wxy);
%vDegree_y_wxy = zeros(1,nEntries_wxy);


% for each pair, perform a deconvolution to obtain w_{x}
if nEntries_hxy > 1
    for i = 1:1:nEntries_hxy - 1
        
        arr_wxy{i} = Deconvolve_Bivariate(arr_hxy{i},arr_hxy{i+1});
        
        % Set the degree of w with respect to x
        %vDegree_x_wxy(i) = vDegree_x_hxy(i) - vDegree_x_hxy(i+1);
        
        % Set the degree of w with respect to y
        %vDegree_y_wxy(i) = vDegree_y_hxy(i) - vDegree_y_hxy(i+1);
        
    end
    % Include the final element of hx in wx
    arr_wxy{i+1} = arr_hxy{i+1};
    
    % Set its degree with respect to x
    %vDegree_x_wxy(i+1) = vDegree_x_hxy(i+1);
    
    % Set its degree with respect to y
    %vDegree_y_wxy(i+1) = vDegree_y_hxy(i+1);
    
    
else
    
    % number of elements in the array h_{i} is equal to one
    arr_wxy{1} = arr_hxy{1};
    
    % Set the degree with respect to x
    %vDegree_x_wxy(1) = vDegree_x_hxy(1);
    
    % Set the degree with respect to y
    %vDegree_y_wxy(1) = vDegree_y_hxy(1);
    
end


end