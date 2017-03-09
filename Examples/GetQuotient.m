function u_root_mult_mat = GetQuotient(f_root_mult_mat, d_root_mult_mat)
% Get the roots of quotient polynomial u(x) given the roots of polynomial 
% f(x), and the roots of polynomial d(x), where d(x) is the GCD of f(x) and
% g(x) and where f(x)/u(x) = d(x)
%
% % Inputs
%
% f_root_mult_mat : (Matrix)
%
% d_root_mult_mat : (Matrix)
%
% % Outputs
%
% u_root_mult_mat : (Matrix)



% Get the number of distinct roots in f(x)
nDistinctRoots_fx = size(f_root_mult_mat, 1);

% Initialise an empty matrix of roots in u(x)
u_root_mult_mat = [];

% Catch the case that the degree of the GCD is zero, and therefore the quotient
% polynomial u(x) is equal to f(x)

[nDistinctRoots_dx] = size(d_root_mult_mat, 1);

if nDistinctRoots_dx == 0
    u_root_mult_mat = f_root_mult_mat;
    return
end

% For each root of f(x), check to see if it exists in d(x).
for i = 1 : 1 : nDistinctRoots_fx
    
    % Get the ith root r_{i}
    root = f_root_mult_mat(i,1);
    
    % Get the multiplicity of the ith root r_{i}
    rootMultiplicity_fx = f_root_mult_mat(i,2);
    
    % Look up the root in roots of d
    if ~isempty(find(d_root_mult_mat(:,1) == root))
    
        % Get the row on which we find the root
        [row_d,~] = find(d_root_mult_mat(:,1) == root);
        
        % Get the multiplicity of the root
        mult_d = d_root_mult_mat(row_d,2);
        
        % Subtract multiplicty in d to obtain multiplicity in quotient
        % polynomial u
        mult_u = rootMultiplicity_fx - mult_d;
        
        % Add the root and its multiplicity to the set of roots for
        % quotient polynomial u
        if mult_u > 0
            u_root_mult_mat = [u_root_mult_mat; root mult_u];
        end
        
    else
        
        u_root_mult_mat = [u_root_mult_mat; root rootMultiplicity_fx];
    end
end