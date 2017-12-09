function [arr_wxy_1, arr_wxy_2, arr_wxy_out] = o_roots_mymethod_xy(arr_wxy_1, arr_wxy_2)
%
% Inputs
%
% arr_wxy_1 : (Array of Matrices)
%
% arr_wxy_2 : (Array of Matrices)
%
% % Outputs
%
% arr_wxy_1 : 
%
% arr_wxy_2 :
%
% arr_wxy_out :


% Get the number of entries in the set of polynomials w(x)
[~,nEntries_wxy_1] = size(arr_wxy_1);

% Get the number of entries in the set of polynomials w(y)
[~,nEntries_wxy_2] = size(arr_wxy_2);


arr_wxy_out = 1;

% For each of the polynomials w(x)
for i = 1:1:nEntries_wxy_1
    
    % Get the degree of w(x) with respect to y.
    [~,m2] = GetDegree_Bivariate(arr_wxy_1{i});
    
    % If the polynomial has a y component (is Bivariate), then must
    % deconvolve with the corresponding polynomial wy_{i}. This will result
    % in a polynomial in x only, a polynomial in y only and the
    % nonseparable polynomial in x and y.
    
    if m2 > 0
        
        [m1,m2] = GetDegree_Bivariate(arr_wxy_1{i});
        [n1,n2] = GetDegree_Bivariate(arr_wxy_2{i});
        
        lowerLimit_t1 = 1;
        upperLimit_t1 = min(m1,n1);
        
        limits_t1 = [lowerLimit_t1, upperLimit_t1];
        
        lowerLimit_t2 = 1;
        upperLimit_t2 = min(m2,n2);
        
        limits_t2 = [lowerLimit_t2, upperLimit_t2];
        
        [fxy_calc_matrix,gxy_calc_matrix, dxy_calc_matrix, uxy_calc_matrix, vxy_calc_matrix,t1,t2]  =...
            o_gcd_mymethod_Bivariate_2Polys(arr_wxy_1{i}, arr_wxy_2{i}, limits_t1, limits_t2);
        
        % Overwrite wx and wy with new values
        % Assign the GCD to the non-Separable part
        arr_wxy_out{i} = dxy_calc_matrix;
        
        % Divide w_{x} by the GCD to obtain w_{x} without y component
        arr_wxy_1{i} = Deconvolve_Bivariate(arr_wxy_1{i},dxy_calc_matrix);
        
        % Divide w_{y} by the GCD to obtain w_{y} without x component
        arr_wxy_2{i} = Deconvolve_Bivariate(arr_wxy_2{i},dxy_calc_matrix);
        
    end
    
end

% %
% For each w_{x,i} get the root
% get the polynomial, whose roots have multiplicty i, in bernstein
% form, where coefficients are in terms of (1-y)^{m-i}y^{i}.


% Given the set of polynomials w(x), get the roots of f(x,y) and their
% corresponding multiplicities.
[root_mult_arr_x] = GetRootMultiplicityArray(arr_wxy_1);
LineBreakLarge()

% Given the set of polynomials w(x), get the roots of f(x,y) and their
% corresponding multiplicities.

% Transpose the entries of w(y)
arr_wxy_2 = cellfun(@transpose,arr_wxy_2,'un',0);

[root_mult_arr_y] = GetRootMultiplicityArray(arr_wxy_2);


display(root_mult_arr_x)
display(root_mult_arr_y)


end


function root_mult_arr_x = GetRootMultiplicityArray(wx)

LineBreakMedium();

fprintf('Roots and Multiplicities of f(x,y) \n')

% Get number of polynomials in the set wx
[~,nEntries_wx] = size(wx);

% Initialise a root and multiplicity matrix;
root_mult_arr_x = [];

for i=1:1:nEntries_wx
    
    % Get the factor
    factor = wx{i};
    
    % Get the number of coefficients
    [r2,~] = size(factor);
    
    if (r2 == 2 )
        aw = wx{i};
        
        % Normalise the factors polynomial coefficients
        aw = aw./aw(1);
        
        % Convert to power form, so that coefficients are in terms of y^{i}
        % rather than (1-y)^{m-i}y^{i}.
        a_pwr = [aw(1,:) ; aw(2,:)-aw(1,:)];
        
        % Obtain the root in terms of y, and set multiplicity to one.
        %fprintf([mfilename ' : ' sprintf('The root of multiplicity %i is given by:\n', i) ]);
        
        
        root_mult_arr_x = [root_mult_arr_x ; -a_pwr(1,:)./a_pwr(2,:) i];
    else
        fprintf([mfilename sprintf('More than one root \n')])
    end
    
end

end