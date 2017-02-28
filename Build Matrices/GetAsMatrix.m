function fxy_matrix = GetAsMatrix(f_vec,m1,m2)
% Given the vector of coefficients of the polynomial f(x,y), 
% format the coefficients as a matrix.

% Initialise an empty matrix fxy
fxy_matrix = zeros(m1+1,m2+1);

% Intialise a counter which will go through each entry of f_vec (The vector
% of coefficients of of f).
count = 1;

% get number of diagonals in the matrix fxy.
num_diags = (m1+1) + (m2+1) -1;

for tot = 0:1:num_diags -1;
    for i = tot:-1:0
        j = tot-i;
        if i > m1 || j> m2
            % restrict to only the i and j values within the matrix.
        else
            fxy_matrix(i+1,j+1) = f_vec(count);
            count = count + 1;
        end
    end
end



end