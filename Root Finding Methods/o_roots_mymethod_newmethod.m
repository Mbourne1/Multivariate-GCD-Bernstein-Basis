function [arr_wxy,vDeg_t_wx] = o_roots_mymethod_newmethod(fxy_matrix,M)
% Given a bivariate polynomial compute the roots by differentiating with
% respect to x.

global SETTINGS

% Set the iteration number
ite = 1;

% Set the first entry of q to be the input polynomial f(x,y)
arr_fxy{1} = fxy_matrix;

% Get the total degree of f(x,y)
vDeg_t_arr_fxy(ite,1) = M(ite);

% Get the degree of f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(arr_fxy{ite});
vDeg_x_arr_fxy(ite,1) = m1;
vDeg_y_arr_fxy(ite,1) = m2;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDeg_x_arr_fxy(ite,1) > 0 && vDeg_y_arr_fxy(ite,1) > 0
    
    
    if (vDeg_x_arr_fxy(ite,1) == 1)
        % The derivative with respect to x is a constant
        
        % The GCD is a constant
        %arr_fxy{ite+1,1} = Differentiate_wrt_x(arr_fxy{ite},vDeg_t_arr_fxy(ite));
        arr_fxy{ite+1,1} = Differentiate_wrt_x(arr_fxy{ite});
       
        % Deconvolve
%         arr_uxy{ite+1,1} = Deconvolve_Bivariate(arr_fxy{ite}, arr_fxy{ite+1},...
%             vDeg_t_arr_fxy(ite), vDeg_t_arr_fxy(ite)-1);
        arr_uxy{ite+1,1} = Deconvolve_Bivariate(arr_fxy{ite}, arr_fxy{ite+1});
        
        
        % Get total degree of d(x,y) and degree with respect to x and y
        vDeg_t_arr_fxy(ite+1,1) = vDeg_t_arr_fxy(ite)-1;
        vDeg_x_arr_fxy(ite+1,1) = vDeg_x_arr_fxy(ite)-1;
        vDeg_y_arr_fxy(ite+1,1) = vDeg_y_arr_fxy(ite);
        
        break;
    end
    
    
    LineBreakLarge()
    fprintf([mfilename ' : ' sprintf('GCD Calculation Loop iteration = %i \n', ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite-1,ite-1)]);
    
    fxy = arr_fxy{ite};
    
    % Get the derivative of f(x,y) with respect to x.
    %gxy = Differentiate_wrt_x(arr_fxy{ite},vDeg_t_arr_fxy(ite));
    gxy = Differentiate_wrt_x(arr_fxy{ite});
    hxy = Differentiate_wrt_y(arr_fxy{ite});
    
    % Get the total degree of f(x,y)
    m =  vDeg_t_arr_fxy(ite);
    
    % Get the total degree of g(x,y)
    n = m - 1;
    
    % Get the total degree of h(x,y)
    o = m - 1;
    
    % Get the upper and lower limit of the degree of the GCD(f, f')
    if ite > 1
        lower_lim = vDeg_t_arr_fxy(ite)-vNumDistinctRoots(ite-1);
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite, lower_lim)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n', ite, upper_lim)]);
    LineBreakLarge();
    
    % GCD is only a scalar with respect to x so set equal to g(x,y).

    limits_t = [lower_lim, upper_lim];
    
    
    [fxy_o, gxy_o, hxy_o, dxy_o, uxy_o, vxy_o, wxy_o, t, t1, t2] = ...
        o_gcd_mymethod_3Polys(fxy, gxy, hxy, m, n, o, limits_t);
    
    %[arr_fxy{ite,1},~,arr_fxy{ite+1,1},arr_uxy{ite,1},arr_vxy{ite,1},t,t1,t2] = o_gcd_mymethod(arr_fxy{ite},gxy,m,n,);
    arr_fxy{ite,1} = fxy_o;
    arr_fxy{ite+1,1} = dxy_o;
    arr_uxy{ite,1} = uxy_o;
    arr_vxy{ite,1} = vxy_o;
    
    
    % Get total degree of d(x,y) and degree with respect to x and y
    
    vDeg_x_arr_fxy(ite+1,1) = t;
    vDeg_y_arr_fxy(ite+1,1) = t1;
    vDeg_t_arr_fxy(ite+1,1) = t2;
    
    % Get number of distinct roots of f(ite)
    vNumDistinctRoots(ite,1) = vDeg_t_arr_fxy(ite) - vDeg_t_arr_fxy(ite+1);
    
    fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDeg_t_arr_fxy(ite+1))]);
    fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,vNumDistinctRoots(ite))]);
    fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite + 1, vDeg_t_arr_fxy(ite+1))])
    
    LineBreakLarge()
    
    % Increment the iteration number
    ite = ite + 1;
    
end

% % Obtain the series h_{i}
% Each h_{x}(i) is obtained by the deconvolution of q_{x}(i) and q_{x}(i+1)

% Get number of polynomials in the series of polynomials q_{i}
[nPolys_fxy] = size(arr_fxy,1);

% %
% %     Get h_{i}(x)
% %

% METHOD refers to method used to compute h{i}(x,y), either by
% deconvolution or from GCD triples computed above.
% From Deconvolutions
% From ux
SETTINGS.HXY_METHOD = 'From Deconvolutions';
switch SETTINGS.HXY_METHOD
    
    case 'From Deconvolutions'
        
        switch SETTINGS.DECONVOLUTION_METHOD
            
            case 'Separate' % Separate deconvolution
                
                for i = 1 : 1 : nPolys_fxy - 1 % For each pair of q_{x}(i) and q_{x}(i+1)
                    
                    % Deconvolve
                    % arr_hxy{i,1} = Deconvolve_Bivariate_Single(arr_fxy{i}, arr_fxy{i+1},vDeg_t_arr_fxy(i), vDeg_t_arr_fxy(i+1));
                    arr_hxy{i,1} = Deconvolve_Bivariate(arr_fxy{i}, arr_fxy{i+1});
                    
                end
            case 'Batch' % Batch deconvolution
                
                % Get the set of polynomials hx{i} from the deconvolution of the
                % set of polynomials fx{i}/fx{i+1}
                
                
                arr_hxy = Deconvolve_Bivariate_Batch(arr_fxy, vDeg_t_arr_fxy);
                
            otherwise
                error([mfilename ' : ' sprintf(' Deconvolution Method is either Separate or Batch')])
                
        end
        
    case 'From ux'
        
        arr_hxy = arr_uxy;
        
    otherwise
        error('err')
        
end

vDeg_x_hxy = vDeg_x_arr_fxy(1:end-1) - vDeg_x_arr_fxy(2:end);
vDeg_y_hxy = vDeg_y_arr_fxy(1:end-1) - vDeg_y_arr_fxy(2:end);
vDeg_t_hxy = vDeg_t_arr_fxy(1:end-1) - vDeg_t_arr_fxy(2:end);


% %
% %
% Obtain the series of polynomials w_{x}{i}

% Each w_{x}(i) is obtained by the deconvolution of h_{x}(i) and h_{x}(i+1)

% Get number of polynomials in the array of h_{x}
[nEntries_hx] = size(arr_hxy,1);

if nEntries_hx > 1
    
    switch SETTINGS.DECONVOLUTION_METHOD
        
        case 'Separate' % Separate deconvolution
            
            for i = 1 : 1 : nEntries_hx - 1 % For each pair of q_{x}(i) and q_{x}(i+1)
                
                % Deconvolve
                %arr_wxy{i,1} = Deconvolve_Bivariate_Single(arr_hxy{i},arr_hxy{i+1},vDeg_t_hxy(i),vDeg_t_hxy(i+1));
                arr_wxy{i,1} = Deconvolve_Bivariate(arr_hxy{i}, arr_hxy{i+1});
            end

        case 'Batch' % Batch deconvolution
            
            arr_wxy = Deconvolve_Bivariate_Batch(arr_hxy, vDeg_t_hxy);
            
        otherwise
            error('err')
            
    end
    
    vDeg_x_wx = vDeg_x_hxy(1:end-1) - vDeg_x_hxy(2:end);
    vDeg_y_wx = vDeg_y_hxy(1:end-1) - vDeg_y_hxy(2:end);
    vDeg_t_wx = vDeg_t_hxy(1:end-1) - vDeg_t_hxy(2:end);
    
    
    % Set the final w_{x}(i+1) to be equal to h_{x}(i+1)
    arr_wxy{end+1,1} = arr_hxy{end};
    
    % Set the final degree structure
    vDeg_x_wx(end) = vDeg_x_hxy(end);
    vDeg_y_wx(end) = vDeg_y_hxy(end);
    vDeg_t_wx(end) = vDeg_t_hxy(end);
    
else
    
    % Number of h_{i} is equal to one, so no deconvolutions to be performed
    arr_wxy{1} = arr_hxy{1};
    
    % Get the degree structure of h_{x,i}
    vDeg_x_wx(1) = vDeg_x_hxy(1);
    vDeg_y_wx(1) = vDeg_y_hxy(1);
    vDeg_t_wx(1) = vDeg_t_hxy(1);
end

for i = 1:1:size(arr_wxy,1)
    fprintf([mfilename ' : ' sprintf('Factors of multiplicity %i :',i) ' \n']);
    factor = arr_wxy{i};
    
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
