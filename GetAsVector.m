function f_vec = GetAsVector(fxy_matrix)
% Given the polynomial f(x,y) in the Bernstein basis, whose coefficients are in
% matrix form, obtain the vector of the coefficients such that the order is
% increasing and the higher power of x is first.

% Get the degree of f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Initialise a counter
count = 1;

% Initialise the vector of coefficients of f(x,y)
f_vec = zeros((m1+1)*(m2+1),1);

% Get number of diagonals in fxy_matrix
nDiags_fxy = (m1+1)+(m2+1)-1;

% For each diagonal of f(x,y), read into the vector, starting lower left to
% upper right.
for tot = 0:1:nDiags_fxy
    
    for i = tot:-1:0
        j = tot - i;
        
        if(i > m1 || j > m2)
            % Do nothing
        else
            
            % Assign matrix entry to next available vector entry
            f_vec(count) = fxy_matrix(i+1,j+1);
            
            % Increment counter
            count = count + 1;
        end
        
    end
    
end


end