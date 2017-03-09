function [max_mtrx,min_mtrx] = GetMaxMin(fxy, n1_k1, n2_k2)
% GetMaxMin(fxy, n1_k1, n2_k2)
%
% Get the maximum and minimum of each entry of f(x,y) in the Sylvester
% matrix.
%
% Inputs
%
% fxy : (Matrix) The Coefficients of polynomial f(x,y) in standard bernstein basis.
%       Given in matrix form so that the rows are in terms of x basis
%       elements and the columns are y basis elements.
%
% n1_k1 : (Int) Degree of v(x,y) with respect to x
%
% n2_k2 : (Int) Degree of v(x,y) with respect to y


% Get the degree of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Take the absolute values of the coefficients
fxy = abs(fxy);

% Build a matrix which stores the maximum values of each coefficient
max_mtrx = zeros(m1+1, m2+1);

% Build a matrix which stores the minimum values of each coefficient
min_mtrx = zeros(m1+1, m2+1);

for i1=0:1:m1
    for i2 = 0:1:m2
        
        % This gives the a_{i1,i2} coefficient of f
        %[maximum,minimum] = GetMaxMin2(fxy(i1+1,i2+1),i1,i2,m1,m2,n1_k1,n2_k2);

        [maximum,minimum] = GetMaxMin2c(fxy(i1+1,i2+1), i1, i2, m1, m2, n1_k1, n2_k2);
        max_mtrx(i1+1,i2+1) = maximum;
        min_mtrx(i1+1,i2+1) = minimum;
        
    end
    
end


end


% % The next three functions compute max and minimum entries by three
% different methods.

function [maximum,minimum] = GetMaxMin2a(a_i1i2,i1,i2,m1,m2,n1_k1,n2_k2)
% 
%
%

% Build a 2 dimensional vector to store all occurences of the coefficient
% a_{i_{1},i_{2}}

val = a_i1i2 .* nchoosek(m1,i1) * nchoosek(m2,i2);

A = ones(n1_k1+1,n2_k2+1) .* val;

v_binoms_numerator_pre_mat = zeros(n1_k1 + 1, 1);
v_binoms_denominator_pre_mat = zeros(n1_k1 + 1, 1);

v_binoms_numerator_post_mat = zeros(n2_k2 + 1, 1);
v_binoms_denominator_post_mat = zeros(n2_k2 + 1, 1);


for j1 = 0:1:n1_k1
    v_binoms_numerator_pre_mat(j1+1) = nchoosek(n1_k1,j1);
    v_binoms_denominator_pre_mat(j1+1) = nchoosek(m1+n1_k1,i1+j1);
end

for j2 = 0:1:n2_k2
    v_binoms_numerator_post_mat(j2+1) = nchoosek(n2_k2,j2);
    v_binoms_denominator_post_mat(j2+1) = nchoosek(m2+n2_k2,i2+j2);
end

global SETTINGS
switch SETTINGS.BOOL_Q
    case 'y'
        pre_mat = diag(v_binoms_numerator_pre_mat./v_binoms_denominator_pre_mat);
        post_mat = diag(v_binoms_numerator_post_mat ./ v_binoms_denominator_post_mat);
    case 'n'
        pre_mat = diag(1./v_binoms_denominator_pre_mat);
        post_mat = diag(1./ v_binoms_denominator_post_mat);
    otherwise
        
end



A = pre_mat * A * post_mat;

% take absolute values of A
A = abs(A);

[max_r,max_c] = find(A==max(A(:)));

[min_r,min_c] = find(A==min(A(:)));


% get the maximum and minimum values. Always use (1) since max or min may
% occur more than once, and we are only interested in one of these values.
maximum = A(max_r(1),max_c(1));
minimum = A(min_r(1),min_c(1));



end


function [maximum,minimum] = GetMaxMin2b(a_i1i2,i1,i2,m1,m2,n1_k1,n2_k2)
% 
%
%

global SETTINGS


% Build a 2 dimensional vector to store all occurences of the coefficient
% a_{i_{1},i_{2}}
A = zeros(n1_k1+1,n2_k2+1);






switch SETTINGS.BOOL_Q 
    case 'y'
        % With the inclusion of Q
        % for each occurence of a_{i_{1},i_{2}} in each column
        for j1 = 0:1:n1_k1
            for j2 = 0:1:n2_k2
                A(j1+1,j2+1) = (a_i1i2 .*  nchoosek(m1,i1) .* nchoosek(m2,i2) .*...
                                nchoosek(n1_k1,j1) .* nchoosek(n2_k2,j2) )./ ...
                                (nchoosek(m1+n1_k1,i1+j1) .* nchoosek(m2+n2_k2,i2+j2));
            end
        end    
    case 'n'
        % Without the inclusion of Q
        % for each occurence of a_{i_{1},i_{2}} in each column
        for j1 = 0:1:n1-k1
            for j2 = 0:1:n2-k2
                A(j1+1,j2+1) = (a_i1i2 .*  nchoosek(m1,i1) .* nchoosek(m2,i2)) ./ ...
                                (nchoosek(m1+n1_k1,i1+j1) .* nchoosek(m2+n2_k2,i2+j2));
            end
        end
end

% take absolute values of A
A = abs(A);

[max_r,max_c] = find(A==max(A(:)));

[min_r,min_c] = find(A==min(A(:)));


% get the maximum and minimum values. Always use (1) since max or min may
% occur more than once, and we are only interested in one of these values.
maximum = A(max_r(1),max_c(1));
minimum = A(min_r(1),min_c(1));



end

function [maximum,minimum] = GetMaxMin2c(a_i1i2,i1,i2,m1,m2,n1_k1,n2_k2)
% 
%
%


% Build a 2 dimensional vector to store all occurences of the coefficient
% a_{i_{1},i_{2}}

val = a_i1i2; 

A = ones(n1_k1+1,n2_k2+1).*val;

for j1 = 0:1:n1_k1
    
    v_binoms_numerator_pre_mat1(j1+1) = nchoosek(i1+j1,i1);
    v_binoms_numerator_pre_mat2(j1+1) = nchoosek(m1+n1_k1-(i1+j1),n1_k1-j1);
    
end

for j2 = 0:1:n2_k2
    
    v_binoms_numerator_post_mat1(j2+1) = nchoosek(i2+j2,i2);
    v_binoms_numerator_post_mat2(j2+1) = nchoosek(m2+n2_k2-(i2+j2),n2_k2-j2);
    
end

binom_denom_pre_mat = nchoosek(m1+n1_k1,n1_k1);
binom_denom_post_mat = nchoosek(m2+n2_k2,n2_k2);

pre_mat = diag(v_binoms_numerator_pre_mat1.*v_binoms_numerator_pre_mat2) ./ binom_denom_pre_mat;
post_mat = diag(v_binoms_numerator_post_mat1.*v_binoms_numerator_post_mat2) ./ binom_denom_post_mat;

A = pre_mat * A * post_mat;

% take absolute values of A
A = abs(A);

[max_r,max_c] = find(A==max(A(:)));

[min_r,min_c] = find(A==min(A(:)));


% get the maximum and minimum values. Always use (1) since max or min may
% occur more than once, and we are only interested in one of these values.
maximum = A(max_r(1),max_c(1));
minimum = A(min_r(1),min_c(1));



end