function fxy_vec = GetAsVector_Symbolic(fxy)
% Given the polynomial f(x,y) in the Bernstein basis, whose coefficients are in
% matrix form, obtain the vector of the coefficients such that the order is
% increasing and the higher power of x is first.
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% % Outputs 
%
% f_vec : (Vector) Vector of ordered symbolic coefficients of f(x,y)


% Get the degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Initialise a counter
count = 1;

% Initialise the vector of coefficients of f(x,y)
fxy_vec = sym(zeros((m1+1)*(m2+1),1));

% Get number of diagonals in fxy_matrix
nDiagonals_fxy = (m1+1)+(m2+1)-1;

% For each diagonal of f(x,y), read into the vector, starting lower left to
% upper right.

method = 'rearranged';

switch method
    case 'standard'
        for tot = 0:1:nDiagonals_fxy
            
            for i = tot:-1:0
                j = tot - i;
                
                if(i > m1 || j > m2)
                    % Do nothing
                else
                    
                    % Assign matrix entry to next available vector entry
                    fxy_vec(count) = fxy(i+1,j+1);
                    
                    % Increment counter
                    count = count + 1;
                end
                
            end
            
        end
    case 'rearranged'
        
        for i = 0:1:m1
            for j = 0:1:m2
                
                fxy_vec(count) = fxy(i+1,j+1);
           
                count = count + 1;
                
            end
        end
        
    otherwise 
        error('err')
end

end