function hxy_matrix = BernMultiply_Bivariate(fxy_matrix,gxy_matrix)
% Given two Bernstein Polynomials where the coefficients are given in a 
% matrix. 
% rows relate to the basis element wrt x 
% cols relate to the basis element wrt y

%               B_{0}(y)  B_{1}(y) ....
%               ________ ________
%    B_{0}(x)  |________|________| ...
%    B_{1}(x)  |________|________| ...
%    B_{2}(x)  |________|________| ...
%       ...       ...       ...

% Get the dimensions and degrees of polynomial fxy
[r, c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get the dimensions and degrees of polynomial gxy
[r, c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

% Get the dimensions and degrees of the product hxy = fxy * gxy
o1 = m1 + n1;
o2 = m2 + n2;


% The coefficient matrix includes the binomial coefficients corresponding
% to a_{i,j} \binom{m_{1}}{i} \binom{m_{2}}{j}
% It is necessary to multiply fxy_matrix by the binomials

% Initialise a vector which multiplies by the rows
bi_m1 = zeros(m1 + 1, 1);
for i = 0 : 1 : m1
    bi_m1(i + 1) = nchoosek(m1, i);
end
bi_m1_mtrx = diag(bi_m1);


% Initialise a vector which multiplies by the columns
bi_m2 = zeros(m2 + 1, 1);
% for each column, multiply by \binom{m2}{i}
for i = 0 : 1 : m2
    bi_m2(i + 1) = nchoosek(m2, i);
end
% form a diagonal matrix
bi_m2_mtrx = diag(bi_m2);

% Since we are multiplying the columns by binomial coefficients 
% \binom{m2}{i} for each column i, we post multiply by the diagonal
% matrix bi_m2.
% Since we are multiplying the rows by binomial coefficients 
% \binom{m1}{i} for each row i, we pre multiply by the diagonal matrix
% bi_m1.

fxy_matrix_bi = bi_m1_mtrx * fxy_matrix * bi_m2_mtrx;


% % Build the matrix T

% Initialise T(f) the matrix whose rows consist of fxy multiplied by the
% basis elements B_{i}(x)B_{j}(y)
nRows_T = (m1 + n1 + 1) * (m2 + n2 + 1);
nCols_T = (n1 + 1)*(n2 + 1);
T = zeros(nRows_T, nCols_T);



% for each diagonal of the matrix g(x,y)
nDiagonals = (n1 + 1) + (n2 + 1) - 1;
count = 1;

for tot = 0:1:nDiagonals
    for i = tot : -1 : 0
        j = tot - i;
    
        % check i is within bounds 0,...,n_{1}
        % check j is within bounds 0,...,
        if i > n1 || j > n2 
        else
            % multiply by the bernstein basis coefficients
            f = mult(fxy_matrix_bi, i, n1, j, n2);
            f_vec = getAsVector(f);
            
            T(:,count) = f_vec;
            count = count + 1;
        end
        
    end
end
    


% % Build the matrix D

count = 1;
nDiagonals = (m1 + n1 + 1) + (m2 + n2 + 1) - 1;
for tot = 0 : 1 : nDiagonals
    for i = tot : -1 : 0
        j = tot - i;
        % check i is within bounds 0,...,n_{1}
        % check j is within bounds 0,...,
        if i > m1 + n1 || j > m2 + n2 
        else
            % 
            d = nchoosek(m1 + n1, i) * nchoosek(m2 + n2, j);
            D(count) = d;
            count = count + 1;
        end
    end
end

D_matrix = diag(1./D);

% % Build the matrix Q

% Initialise the counter
count = 1;

% Get number of diagonals in the matrix gxy_matrix
nDiagonals = (n1 + 1) + (n2 + 1) - 1;

% for each diagonal, read coefficients from left to right, highest power of
% x to minimum power of x.
for tot = 0 : 1 : nDiagonals
    for i = tot : -1 : 0
        j = tot - i;
        if i > n1 || j > n2
        else
            %
            q = nchoosek(n1, i) * nchoosek(n2, j);
            Q(count) = q;
            count = count + 1;
        end
    end
    
end
Q_matrix = diag(Q);

% % Build the coefficient matrix DTQ

DTQ = D_matrix * T * Q_matrix;

% % Build the vector of the coefficients of gxy
g_vec = getAsVector(gxy_matrix);

% % Perform multiplication to obtain hxy in vector form

h_vec = DTQ * g_vec;

% % Convert hxy to matrix form
hxy_matrix = getAsMatrix(h_vec, o1, o2);
end


function f_padd = mult(fxy_matrix, i, n1, j, n2)

[m1, m2] = GetDegree_Bivariate(fxy_matrix);


% get the number of rows and columns in the padded matrix f*basis element
rows = m1 + n1 + 1;
cols = m2 + n2 + 1;

% padd f with zeros
f_padd = zeros(rows,cols);


%fprintf('Multiply %i times with respect to x \n',i)
%fprintf('Multiply %i times with repsect to y \n',j)

i_hat = i + 1;
j_hat = j + 1;

f_padd(i_hat : i + r, j_hat : j + c) = fxy_matrix;


end

function f_vec = getAsVector(fxy_matrix)
% Given the polynomial f in the bernstein basis, whose coefficients are in
% matrix form, obtain the vector of the coefficients such that the order is
% increasing and the higher power of x is first.

[m1,m2] = GetDegree(fxy_matrix);

count = 1;
f_vec = zeros(r*c, 1);

nDiagonals = r + c - 1;


for tot = 0 : 1 : nDiagonals
    for i = tot : -1 : 0
        j = tot - i;
        
        if(i > m1 || j > m2)
        else
            f_vec(count) = fxy_matrix(i + 1, j + 1);
            count = count + 1;
        end
        
    end
    
end


end

function fxy_matrix = getAsMatrix (f_vec, m1, m2)
% Given the vector of coefficients of fxy, format the coefficients as a
% matrix.

% Initialise an empty matrix fxy
fxy_matrix = zeros(m1 + 1, m2 + 1);

count = 1;

% get number of diagonals in the matrix fxy.
nDiagonals = (m1 + 1) + (m2 + 1) -1;

for tot = 0 : 1 : nDiagonals -1
    for i = tot : -1 : 0
        j = tot - i;
        if i > m1 || j> m2
            
        else
            fxy_matrix(i + 1, j + 1) = f_vec(count);
            count = count + 1;
        end
    end
end



end