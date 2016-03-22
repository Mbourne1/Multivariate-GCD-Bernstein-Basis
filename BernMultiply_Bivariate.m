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
[m1,m2] = GetDegree(fxy_matrix);

% Get the dimensions and degrees of polynomial gxy
[n1,n2] = GetDegree(gxy_matrix);

% Get the dimensions and degrees of the product hxy = fxy * gxy
o1 = m1 + n1;
o2 = m2 + n2;


% The coefficient matrix includes the binomial coefficients corresponding
% to a_{i,j} \binom{m_{1}}{i} \binom{m_{2}}{j}
% It is necessary to multiply fxy_matrix by the binomials

% Since we are multiplying the columns by binomial coefficients 
% \binom{m2}{i} for each column i, we post multiply by the diagonal
% matrix bi_m2.
% Since we are multiplying the rows by binomial coefficients 
% \binom{m1}{i} for each row i, we pre multiply by the diagonal matrix
% bi_m1.


fxy_matrix_bi = GetWithBinomials(fxy_matrix);


T = BuildT1(fxy_matrix,n1,n2);

% % Build the matrix D
D = BuildD(0,0,m1,m2,n1,n2);


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

% Get degree strucutre of f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% get the number of rows and columns in the padded matrix f*basis element
rows = m1 + n1 + 1;
cols = m2 + n2 + 1;

% padd f with zeros
f_padd = zeros(rows,cols);


%fprintf('Multiply %i times with respect to x \n',i)
%fprintf('Multiply %i times with repsect to y \n',j)

i_hat = i+1;
j_hat = j+1;

f_padd(i_hat:i+(m1+1),j_hat:j+(m2+1)) = fxy_matrix;


end



