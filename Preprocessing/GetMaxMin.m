function [max_mtrx,min_mtrx] = GetMaxMin(fxy,n1_k1,n2_k2)
% Get the maximum and minimum of each entry of f(x,y) in the Sylvester
% matrix.
%
% Inputs
%
% fxy : The Coefficients of polynomial f(x,y) in standard bernstein basis.
%       Given in matrix form so that the rows are in terms of x basis
%       elements and the columns are y basis elements.
%
% n1_k1 : Degree of v(x,y) with respect to x
%
% n2_k2 : Degree of v(x,y) with respect to y

% Get the degree of polynomial f(x,y)
[m1,m2] = GetDegree(fxy);

% Take the absolute values of the coefficients
fxy = abs(fxy);


% Build a matrix which stores the maximum values of each coefficient
max_mtrx = zeros(m1+1,m2+1);

% Build a matrix which stores the minimum values of each coefficient
min_mtrx = zeros(m1+1,m2+1);

for i1=0:1:m1
    for i2 = 0:1:m2
        % This gives the i1,i2 coefficient of f
        [maximum,minimum] = GetMaxMin2(fxy(i1+1,i2+1),i1,i2,m1,m2,n1_k1,n2_k2);
        
        max_mtrx(i1+1,i2+1) = maximum;
        min_mtrx(i1+1,i2+1) = minimum;
        
    end
    
end


end

function [maximum,minimum] = GetMaxMin2(a_i1i2,i1,i2,m1,m2,n1_k1,n2_k2)
% Note this function assumes the inclusion of Q in the coefficient matrix.

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