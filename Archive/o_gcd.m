function [] = o_gcd()

% Performing a polynomial multiplication


% Get a polynomial f

% Let m1 be the degree of polynomial f in the first bernstein basis
m1 = 2;

% Let m2 be the degree of polynomial f in the second bernstein basis
m2 = 2;

% Let n1 be the degree of polynomial g in the first bernstein basis
% (x)(1-x)
n1 = 3;

% Let n2 be the degree of polynomial g in the second bernstein basis y
% (y-1)
n2 = 3;

% let the coefficients of f be given by
% Build matrix of coefficients of polynomial f where

matrix_f = sym('a',[m2+1 m1+1])

for i = 0:1:m1
    for j = 0:1:m2
        matrix_f(i+1,j+1) = matrix_f(i+1,j+1)* sym(nchoosek(2,i)) * sym(nchoosek(2,j));
    end
end

f_vec = [];
for i = -10 : 1 : +10
    f_vec = [f_vec ; flipud(diag(flipud(matrix_f),i))];
end
f_vec;

% Build matrix of coefficients of polynomial G

matrix_g = sym('b',[n2+1 n1+1]);

for i = 0:1:n1
    for j = 0:1:n2
        matrix_g(i+1,j+1) = matrix_g(i+1,j+1)* sym(nchoosek(3,i)) * sym(nchoosek(3,j));
    end
end

matrix_g;


% Build an empty matrix up to the highest powers m1+n1 and m2+n2

g_matrix = sym(zeros(m1+n1+1,m2+n2+1))


% Get number of multiplications of g
(m1+1) * (m2+2);


% get g* x^0 (1-x)^2 * y^0(1-y)^2

for i = 0:1:n1
    for j = 0:1:n2
        g_matrix(i+1,j+1) = matrix_g(i+1,j+1);
    end
end

g_matrix;

g_matrix_1 = g_matrix
g_matrix_2 = circshift(g_matrix,[0,1]);
g_matrix_3 = circshift(g_matrix,[1,0]);
g_matrix_4 = circshift(g_matrix,[0,2]);
g_matrix_5 = circshift(g_matrix,[1,1]);
g_matrix_6 = circshift(g_matrix,[2,0]);
g_matrix_7 = circshift(g_matrix,[1,2]);
g_matrix_8 = circshift(g_matrix,[2,1]);
g_matrix_9 = circshift(g_matrix,[2,2]);

% Convert these matrices to column vectors
g_vec_1 = [];
g_vec_2 = [];
g_vec_3 = [];
g_vec_4 = [];
g_vec_5 = [];
g_vec_6 = [];
g_vec_7 = [];
g_vec_8 = [];
g_vec_9 = [];


g_vec_1 = [];
for i = -10 : 1 : +10
    g_vec_1 = [g_vec_1 ; flipud(diag(flipud(g_matrix_1),i))];
    g_vec_2 = [g_vec_2 ; flipud(diag(flipud(g_matrix_2),i))];
    g_vec_3 = [g_vec_3 ; flipud(diag(flipud(g_matrix_3),i))];
    g_vec_4 = [g_vec_4 ; flipud(diag(flipud(g_matrix_4),i))];
    g_vec_5 = [g_vec_5 ; flipud(diag(flipud(g_matrix_5),i))];
    g_vec_6 = [g_vec_6 ; flipud(diag(flipud(g_matrix_6),i))];
    g_vec_7 = [g_vec_7 ; flipud(diag(flipud(g_matrix_7),i))];
    g_vec_8 = [g_vec_8 ; flipud(diag(flipud(g_matrix_8),i))];
    g_vec_9 = [g_vec_9 ; flipud(diag(flipud(g_matrix_9),i))];
    
end

matrix = [...
    g_vec_1...
    g_vec_2...
    g_vec_3...
    g_vec_4...
    g_vec_5...
    g_vec_6...
    g_vec_7...
    g_vec_8...
    g_vec_9...
    ]


matrix * f_vec



end

