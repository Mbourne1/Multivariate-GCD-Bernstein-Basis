function [wx,wy,wxy] = o_roots_mymethod_xy(wx,wy,vDegt_wx,vDegt_wy)
%
% Inputs
%
% wx : set of polynomials w_{i} where w_{i} is a factor of multiplicity i
% in f(x,y).
%
%

% Get the number of entries in the set of polynomials w(x)
[~,num_entries_wx] = size(wx);

% Get the number of entries in the set of polynomials w(y)
[~,num_entries_wy] = size(wy);


wxy = 1;

% For each of the polynomials w(x)
for i = 1:1:num_entries_wx
    
    % Get the degree of w(x) with respect to y.
    [~,m2] = GetDegree(wx{i});
    
    % If the polynomial has a y component (is Bivariate), then must
    % deconvolve with the corresponding polynomial wy_{i}. This will result
    % in a polynomial in x only, a polynomial in y only and the
    % nonseparable polynomial in x and y.
    
    if m2 > 0
        
        lower_lim = 1;
        upper_lim = min(vDegt_wx(i),vDegt_wy(i));
        
        [fxy_calc_matrix,gxy_calc_matrix, dxy_calc_matrix, uxy_calc_matrix, vxy_calc_matrix,t,t1,t2]  =...
            o_gcd_mymethod(wx{i},wy{i},vDegt_wx(i),vDegt_wy(i),[lower_lim, upper_lim]);
        
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


% Given the set of polynomials w(x), get the roots of f(x,y) and their
% corresponding multiplicities.
[root_mult_arr_x] = GetRootMultiplicityArray(wx);
LineBreakLarge()

% Given the set of polynomials w(x), get the roots of f(x,y) and their
% corresponding multiplicities.

% Transpose the entries of w(y)
wy = cellfun(@transpose,wy,'un',0);

[root_mult_arr_y] = GetRootMultiplicityArray(wy);


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