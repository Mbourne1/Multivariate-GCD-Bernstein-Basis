function hxy_matrix = BernMultiply_Bivariate(fxy_matrix,gxy_matrix)
% Given two Bernstein Polynomials where the coefficients are given in a 
% matrix 
% rows relate to the basis element wrt x 
% cols relate to the basis element wrt y

%               B_{0}(y)  B_{1}(y) ....
%               ________ ________
%    B_{0}(x)  |________|________| ...
%    B_{1}(x)  |________|________| ...
%    B_{2}(x)  |________|________| ...
%       ...       ...       ...

% Get the dimensions and degrees of polynomial fxy
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get the dimensions and degrees of polynomial gxy
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

% Get the dimensions and degrees of the product hxy = fxy * gxy
o1 = m1 + n1;
o2 = m2 + n2;


% The coefficient matrix includes the binomial coefficients corresponding
% to a_{i,j} \binom{m_{1}}{i} \binom{m_{2}}{j}
% It is necessary to multiply fxy_matrix by the binomials

% Initialise a vector which multiplies by the rows
bi_m1 = GetBinomials(m1);
bi_m1_mtrx = diag(bi_m1);


% Initialise a vector which multiplies by the columns
bi_m2 = GetBinomials(m2);
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
num_rows_T = (m1+n1+1)*(m2+n2+1);
num_cols_T = (n1+1)*(n2+1);
T = zeros(num_rows_T,num_cols_T);



% for each diagonal of the matrix g(x,y)
num_diags = (n1+1) + (n2+1) -1;
count = 1;
for tot = 0:1:num_diags
    for i = tot:-1:0
        j = tot - i;
    
        % check i is within bounds 0,...,n_{1}
        % check j is within bounds 0,...,
        if i > n1 || j > n2 
        else
            % multiply by the bernstein basis coefficients
            f = mult(fxy_matrix_bi,i,n1,j,n2);
            f_vec = GetAsVector(f);
            T(:,count) = f_vec;
            count = count + 1;
        end
        
    end
end
    


% % Build the matrix D
Pre_matrix = diag(GetBinomials(m1+n1));
Post_matrix =  diag(GetBinomials(m2+n2));
D = Pre_matrix * ones(m1+n1+1,m2+n2+1) * Post_matrix;
D = diag(GetAsVector(1./D));

% % Build the matrix Q
Pre_matrix = diag(GetBinomials(n1));
Post_matrix = diag(GetBinomials(n2));
temp =  Pre_matrix * ones(n1+1,n2+1) * Post_matrix;
Q = diag(GetAsVector(temp));



% % Build the coefficient matrix DTQ

DTQ = D * T * Q;

% % Build the vector of the coefficients of gxy
g_vec = GetAsVector(gxy_matrix);

% % Perform multiplication to obtain hxy in vector form

h_vec = DTQ * g_vec;

% % Convert hxy to matrix form
hxy_matrix = GetAsMatrix(h_vec,o1,o2);
end


function f_padd = mult(fxy_matrix,i,n1,j,n2)

[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% get the number of rows and columns in the padded matrix f*basis element
rows = m1 + n1 + 1;
cols = m2 + n2 + 1;

% padd f with zeros
f_padd = zeros(rows,cols);


%fprintf('Multiply %i times with respect to x \n',i)
%fprintf('Multiply %i times with repsect to y \n',j)

i_hat = i+1;
j_hat = j+1;

f_padd(i_hat:i+r,j_hat:j+c) = fxy_matrix;


end



