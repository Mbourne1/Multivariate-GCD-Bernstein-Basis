function d_roots = getDivisor(f_roots,g_roots)
% Given the set of roots of f(x) and roots of f(x) and their multiplicities,
% get the common roots of polynomials f(x) and g(x) given by d(x).

% get the number of roots in polynomial f
num_roots_f = size(f_roots,1);

% Initialise the set of roots of d(x)
d_roots = [];

% for each root in f(x), check to see if it exists in g(x)
for i = 1:1:num_roots_f
    
    % Get the root of f(x)
    root = f_roots(i,1);
    
    % Get the multiplicity of root r_{i}
    mult_root_in_f = f_roots(i,2);
    
    % Get the number of distinct roots in g
    [distinct_roots_g,~] = size(g_roots);
    if  distinct_roots_g == 0
        return
    end
    
    % Look if the root r_{i} exists in g(x)
    if ~isempty(find(g_roots(:,1) == root));
        
        % Get the index of the row which corresponds to the root r_{i} in
        % the matrix of roots of g(x)
        [row_d,~] = find(g_roots(:,1) == root);
        
        % Get the multiplicity of the root r_{i} in g(x)
        mult_root_in_g = g_roots(row_d,2);
        
        % Calculate the multiplicity of the root in d(x)
        mult_root_in_d = min(mult_root_in_f,mult_root_in_g);
        
        % Add the root to d(x)
        d_roots = [d_roots ; root mult_root_in_d]; 
    end
end


end

