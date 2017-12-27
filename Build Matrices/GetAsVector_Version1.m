function f_vec = GetAsVector_Version1(fxy)

% Get the degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Initialise a counter
count = 1;

% Initialise the vector of coefficients of f(x,y)
f_vec = zeros((m1 + 1)*(m2 + 1),1);

% Get number of diagonals in fxy_matrix
nDiagonals_fxy = (m1 + 1)+(m2 + 1)-1;

% For each diagonal of f(x,y), read into the vector, starting lower left to
% upper right.
for tot = 0 : 1 : nDiagonals_fxy
    
    for i1 = tot : -1 : 0
        i2 = tot - i1;
        
        if(i1 > m1 || i2 > m2)
            % Do nothing
        else
            
            % Assign matrix entry to next available vector entry
            f_vec(count) = fxy(i1+1, i2+1);
            
            % Increment counter
            count = count + 1;
        end
        
    end
    
end