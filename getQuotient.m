function u_roots = getQuotient(f_roots,d_roots)
% Fet the roots of quotient polynomial u(x) given the roots of polynomial 
% f(x), and the roots of polynomial d(x), where d(x) is the GCD of f(x) and
% g(x) and where f(x)/u(x) = d(x)

% Get the number of distinct roots in f(x)
nRoots_f_x = size(f_roots,1);

% Initialise an empty matrix of roots in u(x)
u_roots = [];

% Catch the case that the degree of the GCD is zero, and therefore the quotient
% polynomial u(x) is equal to f(x)
[r,~] = size(d_roots);
if r == 0
    u_roots = f_roots;
    return
end

% For each root of f(x), check to see if it exists in d(x).
for i = 1:1:nRoots_f_x
    
    % Get the root r_{i}
    root = f_roots(i,1);
    
    % Get the multiplicity of the root r_{i}
    mult_f = f_roots(i,2);
    
    % Look up the root in roots of d
    if ~isempty(find(d_roots(:,1) == root));
    
        % Get the row on which we find the root
        [row_d,~] = find(d_roots(:,1) == root);
        
        % Get the multiplicity of the root
        mult_d = d_roots(row_d,2);
        
        % Subtract multiplicty in d to obtain multiplicity in quotient
        % polynomial u
        mult_u = mult_f - mult_d;
        
        % Add the root and its multiplicity to the set of roots for
        % quotient polynomial u
        if mult_u > 0
            u_roots = [u_roots; root mult_u];
        end
        
    else
        
        u_roots = [u_roots; root mult_f];
    end
end