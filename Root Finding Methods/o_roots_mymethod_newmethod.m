function [arr_wxy] = o_roots_mymethod_newmethod(fxy_matrix)
% Given a bivariate polynomial compute the roots by differentiating with
% respect to x.
%
% % Inputs
%
% fxy_matrix : (Matrix) Coefficients of the polynomial f(x,y)
%
% % Outputs
%
% arr_wxy : (Array of Matrices) Array containing the coefficients of the
% set of polynomials w_{i}(x,y) where each w_{i}(x,y) is the product of
% factors of f(x,y) with multiplicity i.


% Get the array of polynomials f_{i}(x,y)
arr_fxy = GetArray_fxy(fxy_matrix, 'Total');

% Get the array of polynomials h_{i}(x,y)
arr_hxy = GetArray_hxy(arr_fxy);

% Get the array of polynomials w_{i}(x,y)
arr_wxy = GetArray_wxy(arr_hxy);


% Obtain the series of polynomials w_{x}{i}

% Get number of polynomials in array w_{i}(x,y)
nPolys_wxy = length(arr_wxy);

for i = 1 : 1 : nPolys_wxy
    
    fprintf([mfilename ' : ' sprintf('Factors of multiplicity %i :',i) ' \n']);
    
    % Get the ith factor
    factor = arr_wxy{i};
    
    % 
    if (length(factor) > 1)
        
        
        fprintf('\n')
        disp(factor./factor(2));
        %disp(factor);
        
    else
        
        fprintf('\n')
        display(1)
        
    end
    
    
end

end



function arr_fxy = GetArray_fxy(fxy_matrix, degree_computation_method)
%
% % Inputs
%
% fxy_matrix : (Matrix) Coefficients of the polynomial f(x,y)
%
% degree_computation_method : (String) 
%   All Subresutlants
%   Linear Method : 
%   
%   
%
% % Outputs
%
% arr_fxy : (Array of Matrices) Matrices containing coefficients of the set
% of polynomials f_{i}(x,y)


% Set the iteration number
ite = 1;

% Set the first entry of q to be the input polynomial f(x,y)
arr_fxy{1} = fxy_matrix;

% Get the degree of f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(arr_fxy{ite});

vDegree_x_arr_fxy(ite,1) = m1;
vDegree_y_arr_fxy(ite,1) = m2;


% Initial rank range
rank_range = [-16 0];

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDegree_x_arr_fxy(ite,1) > 0 && vDegree_y_arr_fxy(ite,1) > 0
    
    
    if (vDegree_x_arr_fxy(ite,1) == 1)
        % The derivative with respect to x is a constant
        
        % The GCD is a constant
        arr_fxy{ite+1, 1} = Differentiate_wrt_x(arr_fxy{ite});
       
        % Deconvolve
        arr_uxy{ite+1, 1} = Deconvolve_Bivariate(arr_fxy{ite}, arr_fxy{ite+1});
        
        
        % Get degree structure of d(x,y)
        vDegree_x_arr_fxy(ite+1, 1) = vDegree_x_arr_fxy(ite)-1;
        vDegree_y_arr_fxy(ite+1, 1) = vDegree_y_arr_fxy(ite);
         
        break;
    end
    
    
    LineBreakLarge()
    fprintf([mfilename ' : ' sprintf('GCD Calculation Loop iteration = %i \n', ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite-1,ite-1)]);
    
    fxy = arr_fxy{ite};
    
    % Get the derivative of f(x,y) with respect to x.
    gxy = Differentiate_wrt_x(arr_fxy{ite});
    hxy = Differentiate_wrt_y(arr_fxy{ite});
    
    
    % Get the relative degree of f(x,y) g(x,y) and h(x,y)
    [m1, m2] = GetDegree_Bivariate(fxy);
    [n1, n2] = GetDegree_Bivariate(gxy);
    [o1, o2] = GetDegree_Bivariate(hxy);
    
    % Get the upper and lower limit of the degree of the GCD(f, f')
    if ite > 1
        
        lowerLimit_t1 = vDegree_x_arr_fxy(ite) - vNumberDistinctRoots_x(ite-1);
        upperLimit_t1 = min([m1, n1, o1]);
        
        lowerLimit_t2 = vDegree_y_arr_fxy(ite) - vNumberDistinctRoots_y(ite-1);
        upperLimit_t2 = min([m2, n2, o2]);
        
    else

        
        lowerLimit_t1 = 0;
        upperLimit_t1 = min([m1, n1, o1]);
        
        lowerLimit_t2 = 0;
        upperLimit_t2 = min([m2, n2, o2]);
    end
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i} with respect to x: %i \n', ite, lowerLimit_t1)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i} with respect to x: %i \n', ite, upperLimit_t1)]);
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i} with respect to y: %i \n', ite, lowerLimit_t2)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i} with respect to y: %i \n', ite, upperLimit_t2)]);
    
    LineBreakLarge();
    
    % GCD is only a scalar with respect to x so set equal to g(x,y).
    limits_t1 = [lowerLimit_t1, upperLimit_t1];
    limits_t2 = [lowerLimit_t2, upperLimit_t2];
    
    
    
    
    % Compute GCD
    [fxy_o, gxy_o, hxy_o, dxy_o, uxy_o, vxy_o, wxy_o, t1, t2, rank_range] = ...
        o_gcd_mymethod_Bivariate_3Polys(fxy, gxy, hxy, limits_t1, limits_t2, rank_range, degree_computation_method);
    
    
    
    %[arr_fxy{ite,1},~,arr_fxy{ite+1,1},arr_uxy{ite,1},arr_vxy{ite,1},t,t1,t2] = o_gcd_mymethod(arr_fxy{ite},gxy,m,n,);
    arr_fxy{ite, 1} = fxy_o;
    arr_fxy{ite + 1, 1} = dxy_o;
    
    % Get total structure of d(x,y)
    vDegree_x_arr_fxy(ite+1, 1) = t1;
    vDegree_y_arr_fxy(ite+1, 1) = t2;
    
    % Get number of distinct roots of f(ite)
    vNumberDistinctRoots_x(ite, 1) = vDegree_x_arr_fxy(ite) - vDegree_x_arr_fxy(ite+1);
    vNumberDistinctRoots_y(ite, 1) = vDegree_y_arr_fxy(ite) - vDegree_y_arr_fxy(ite+1);
    
    LineBreakLarge()
    fprintf([mfilename ' : ' sprintf('The degree of the GCD with respect to x is : %i \n',  vDegree_x_arr_fxy(ite+1,1))]);
    fprintf([mfilename ' : ' sprintf('The degree of the GCD with respect to y is : %i \n',  vDegree_y_arr_fxy(ite+1,1))]);
    LineBreakLarge()
    
    % Increment the iteration number
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



% Get number of polynomials in array f_{i}(x,y)
nPolys_fxy = length(arr_fxy);
nPolys_hxy = nPolys_fxy - 1;

global SETTINGS
SETTINGS.HXY_METHOD = 'From Deconvolutions';

% Initialise array to store polynomials h_{i}(x,y)
arr_hxy = cell(nPolys_hxy);



switch SETTINGS.DECONVOLUTION_METHOD

    case 'Separate' % Separate deconvolution

        for i = 1 : 1 : nPolys_hxy

            arr_hxy{i,1} = Deconvolve_Bivariate(arr_fxy{i}, arr_fxy{i+1});

        end

    case 'Batch' % Batch deconvolution

        % Get the set of polynomials hx{i} from the deconvolution of the
        % set of polynomials fx{i}/fx{i+1}

        error('Code not complete');
        arr_hxy = Deconvolve_Bivariate_Batch(arr_fxy, vDeg_t_arr_fxy);

    otherwise
        error([mfilename ' : ' sprintf(' Deconvolution Method is either Separate or Batch')])

end
        
 

%vDegree_x_hxy = vDegree_x_arr_fxy(1 : end - 1) - vDegree_x_arr_fxy(2 : end);
%vDegree_y_hxy = vDegree_y_arr_fxy(1 : end - 1) - vDegree_y_arr_fxy(2 : end);



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



% 
% Each w_{x}(i) is obtained by the deconvolution of h_{x}(i) and h_{x}(i+1)

% Get number of polynomials in the array of h_{x}
[nPolysArray_hxy] = size(arr_hxy, 1);


global SETTINGS

if nPolysArray_hxy > 1
    
    switch SETTINGS.DECONVOLUTION_METHOD
        
        case 'Separate' % Separate deconvolution
            
            for i = 1 : 1 : nPolysArray_hxy - 1 % For each pair of q_{x}(i) and q_{x}(i+1)
                
                % Deconvolve
               
                arr_wxy{i,1} = Deconvolve_Bivariate(arr_hxy{i}, arr_hxy{i+1});
            end

        case 'Batch' % Batch deconvolution
            
            error('Code not complete')
            arr_wxy = Deconvolve_Bivariate_Batch(arr_hxy);
            
        otherwise
            error('err')
            
    end
    
    %vDegree_x_wxy = vDegree_x_hxy(1:end-1) - vDegree_x_hxy(2:end);
    %vDegree_y_wxy = vDegree_y_hxy(1:end-1) - vDegree_y_hxy(2:end);
  
    
    % Set the final w_{x}(i+1) to be equal to h_{x}(i+1)
    arr_wxy{end+1,1} = arr_hxy{end};
    
    % Set the final degree structure
    %vDegree_x_wxy(end) = vDegree_x_hxy(end);
    %vDegree_y_wxy(end) = vDegree_y_hxy(end);
  
    
else
    
    % Number of h_{i} is equal to one, so no deconvolutions to be performed
    arr_wxy{1} = arr_hxy{1};
    
    % Get the degree structure of h_{x,i}
    %vDegree_x_wxy(1) = vDegree_x_hxy(1);
    %vDegree_y_wxy(1) = vDegree_y_hxy(1);
  
end

end
