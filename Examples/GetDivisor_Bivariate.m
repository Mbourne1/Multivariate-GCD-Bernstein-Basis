function d_root_mult_mat = GetDivisor_Bivariate(f_root_mult_mat, g_root_mult_mat)
% Given the set of roots of f(x) and roots of f(x) and their multiplicities,
% get the common roots of polynomials f(x) and g(x) given by d(x).
%
% f_root_mult : 
%


% Get the number of distinct factors in polynomial f(x)
nDistinctFactors_fx = size(f_root_mult_mat,1);

% Initialise the set of roots of d(x)
d_root_mult_mat = [];

% For each distinct factor in f(x), check to see if it exists in g(x)
for i = 1 : 1 : nDistinctFactors_fx
    
    % Get the root of f(x)
    root = f_root_mult_mat(i, 1);
    
    % Get the multiplicity of root r_{i}
    rootMultiplicity_in_fx = f_root_mult_mat(i, 2);
    
    % Get the number of distinct roots in g
    [nDistinct_roots_g] = size(g_root_mult_mat, 1);
    
    if  nDistinct_roots_g == 0
        return
    end
    
    % Look if the root r_{i} exists in g(x)
    if ~isempty(find(g_root_mult_mat(:,1) == root))
        
        % Get the index of the row which corresponds to the root r_{i} in
        % the matrix of roots of g(x)
        [row_d,~] = find(g_root_mult_mat(:,1) == root);
        
        % Get the multiplicity of the root r_{i} in g(x)
        mult_root_in_g = g_root_mult_mat(row_d,2);
        
        % Calculate the multiplicity of the root in d(x)
        mult_root_in_d = min(rootMultiplicity_in_fx,mult_root_in_g);
        
        % Add the root to d(x)
        d_root_mult_mat = [d_root_mult_mat ; root mult_root_in_d]; 
    end
end


end

