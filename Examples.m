function [fxy_matrix_exact,gxy_matrix_exact,dxy_matrix_exact,...
    uxy_matrix_exact,vxy_matrix_exact,...
    m,m1,m2,...
    n,n1,n2,...
    t,t1,t2] = Examples(ex_num)

switch ex_num
    
    case '-1'
        f_x_roots = [
            2   1;
            3   1;
            ];
        f_y_roots = [
            2   1;
            5   1;
            ];
        
        g_x_roots = [
            2   1;
            1   1;
            ];
        g_y_roots = [
            2   1;
            6   1;
            ];
        
        d_x_roots = [
            2   1;
            ];
        
        d_y_roots = [
            2   1;
            ];
        
        
    case '0'
        f_x_roots = [
            0.2   1;
            0.3   1;
            ];
        f_y_roots = [
            0.2   1;
            0.3   1;
            ];
        
        g_x_roots = [
            0.2   1;
            0.7   1;
            ];
        g_y_roots = [
            0.2   1;
            0.7   1;
            ];
        
        d_x_roots = [
            0.2   1;
            ];
        
        d_y_roots = [
            0.2   1;
            ];
        
    case '1'
        f_x_roots = [
            0.2   1;
            0.3   1;
            0.4   1;
            0.9   1;
            ];
        f_y_roots = [
            0.2   1;
            0.3   1;
            0.4   1;
            0.5   1;
            0.6   1;
            ];
        
        g_x_roots = [
            0.2   1;
            0.4   1;
            0.7   1;
            0.9   1;
            ];
        g_y_roots = [
            0.2   1;
            0.5   1;
            0.8   1;
            ];
        
        d_x_roots = [
            0.2   1;
            0.4   1;
            0.9   1;
            ];
        
        d_y_roots = [
            0.2   1;
            0.5   1;
            ];
        
        
        
        
    case '3'
        f_x_roots = [
            0.3   1;
            0.4   1;
            0.9   1;
            ];
        f_y_roots = [
            0.2   1;
            0.3   1;
            0.5   1;
            ];
        
        g_x_roots = [
            0.3   1;
            0.4   1;
            ];
        
        g_y_roots = [
            0.1   1;
            0.5   1;
            0.8   1;
            ];
        
        d_x_roots = [
            0.3   1;
            0.4   1;
            ];
        
        d_y_roots = [
            0.5   1;
            ];
        
    case '4'
        f_x_roots = [
            0.3   2;
            0.4   3;
            0.9   1;
            ];
        f_y_roots = [
            0.2   1;
            0.3   1;
            0.5   1;
            ];
        
        g_x_roots = [
            0.3   2;
            0.4   3;
            ];
        
        g_y_roots = [
            0.1   1;
            0.5   1;
            0.8   1;
            ];
        
        d_x_roots = [
            0.3   2;
            0.4   3;
            ];
        
        d_y_roots = [
            0.5   1;
            ];
        
    case '5'
        f_x_roots = [
            0.3   2;
            0.4   3;
            0.9   1;
            0.1   1;
            ];
               
        f_y_roots = [
            0.2   1;
            0.3   1;
            0.5   1;
            ];
        
        g_x_roots = [
            0.3   2;
            0.4   3;
            0.1   1;
            ];
        
        g_y_roots = [
            0.1   1;
            0.5   1;
            0.8   1;
            ];
        
        d_x_roots = [
            0.3   2;
            0.4   3;
            0.1   1;
            ];
        
        d_y_roots = [
            0.5   1;
            ];
        
    case '6'
        f_x_roots = [
            0.3     2;
            0.4     3;
            0.9     1;
            0.1     1;
            1.1     1;
            ];
        f_y_roots = [
            0.2   1;
            0.3   1;
            0.5   1;
            ];
        
        g_x_roots = [
            0.3   2;
            0.4   3;
            0.1   1;
            ];
        
        g_y_roots = [
            0.1   1;
            0.5   1;
            0.8   1;
            ];
        
        d_x_roots = [
            0.3   2;
            0.4   3;
            0.1   1;
            ];
        
        d_y_roots = [
            0.5   1;
            ];
        
    case '7'
        f_x_roots = [
            0.3     2;
            0.4     3;
            0.9     1;
            0.1     1;
            1.1     1;
            1.51296     1;
            ];
        f_y_roots = [
            0.2   1;
            0.3   1;
            0.5   1;
            ];
        
        g_x_roots = [
            0.3   2;
            0.4   3;
            0.1   1;
            ];
        
        g_y_roots = [
            0.1   1;
            0.5   1;
            0.8   1;
            ];
        
        d_x_roots = [
            0.3   2;
            0.4   3;
            0.1   1;
            ];
        
        d_y_roots = [
            0.5   1;
            ];
    case '99'
        f_x_roots = [
            2   1;
            3   1;
            4   1;
            9   1;
            ];
        f_y_roots = [
            2   1;
            3   1;
            4   1;
            5   1;
            6   1;
            ];
        
        g_x_roots = [
            2   1;
            4   1;
            7   1;
            9   1;
            ];
        g_y_roots = [
            2   1;
            5   1;
            7   1;
            ];
        
        d_x_roots = [
            2   1;
            4   1;
            9   1;
            ];
        
        d_y_roots = [
            2   1;
            5   1;
            ];
end

u_x_roots = getQuotient(f_x_roots,d_x_roots);
v_x_roots = getQuotient(g_x_roots,d_x_roots);
u_y_roots = getQuotient(f_y_roots,d_y_roots);
v_y_roots = getQuotient(g_y_roots,d_y_roots);

%% Build polynomial d(x,y)
% Build the Polynomial of the GCD in x terms and y terms
d_x_poly = BuildPoly(d_x_roots);
d_y_poly = BuildPoly(d_y_roots);

% Get matrix of polynomial coefficients including the binomial coefficients
dxy_matrix_bi = d_x_poly' * d_y_poly;

% Get degree of d in terms of x
d1 = length(d_x_poly)-1;
% Get degree of d in terms of y
d2 = length(d_y_poly)-1;

% Print out information about the polynomial d.
fprintf('Input Polynomial d: \n')
fprintf('Degree with respect to x : \n')
fprintf('  d1 = %i \n',d1)
fprintf('Degree with respect to y : \n')
fprintf('  d2 = %i \n',d2)
fprintf('\n')

% Strip the binomial coefficients from d_bi in the scaled bernstein basis,
% to obtain d in the standard Bernstein Basis
dxy_matrix_exact = zeros(d1+1,d2+1);
for i = 0:1:d1
    dxy_matrix_exact(i+1,:) = dxy_matrix_bi(i+1,:) ./ nchoosek(d1,i);
end
for j = 0:1:d2
    dxy_matrix_exact(:,j+1) = dxy_matrix_exact(:,j+1) ./ nchoosek(d2,j);
end


%%      Build polynomial f(x,y)

% Build Polynomial f in x terms and y terms
f_x_poly_bi = BuildPoly(f_x_roots);
f_y_poly_bi = BuildPoly(f_y_roots);
% Build the polynomial f(x,y)
fxy_matrix_bi =  f_x_poly_bi' * f_y_poly_bi;

% Get degree of f in terms of x
m1 = length(f_x_poly_bi)-1;
% Get degree of f in terms of y
m2 = length(f_y_poly_bi)-1;

fprintf('Input Polynomial f: \n')
fprintf('Degree with respect to x : \n')
fprintf('  m1 = %i \n',m1)
fprintf('Degree with respect to y : \n')
fprintf('  m2 = %i \n',m2)
fprintf('\n')

% Strip the binomial coefficients from f_bi in the scaled Bernstein basis
% to obtain the coefficients of f in the standard Bernstein Basis)
fxy_matrix = zeros(m1+1,m2+1);
for i = 0:1:m1
    fxy_matrix(i+1,:) = fxy_matrix_bi(i+1,:) ./ nchoosek(m1,i);
end
for j = 0:1:m2
    fxy_matrix(:,j+1) = fxy_matrix(:,j+1) ./ nchoosek(m2,j);
end

%%          Build Polynomial g(x,y)

% Build Polynomial g in x terms and y terms
g_x_poly_bi = BuildPoly(g_x_roots);
g_y_poly_bi = BuildPoly(g_y_roots);

% Get matrix of polynomial g(x,y) coefficients including the binomial
% coefficients.
gxy_matrix_bi =  g_x_poly_bi' * g_y_poly_bi;

% Get degree of f in terms of x
n1  = length(g_x_poly_bi)-1;
% Get degree of f in terms of y
n2 = length(g_y_poly_bi)-1;

fprintf('Input Polynomial g: \n')
fprintf('Degree with respect to x : \n')
fprintf('  n1 = %i \n',n1)
fprintf('Degree with respect to y : \n')
fprintf('  n2 = %i \n',n2)
fprintf('\n')

% Strip the binomial coefficients from g_bi in the scaled Bernstein Basis,
% to obtain g in the standard Bernstein Basis
gxy_matrix = zeros(n1+1,n2+1);
for i = 0:1:n1
    gxy_matrix(i+1,:) = gxy_matrix_bi(i+1,:) ./ nchoosek(n1,i);
end
for j = 0:1:n2
    gxy_matrix(:,j+1) = gxy_matrix(:,j+1) ./ nchoosek(n2,j);
end

%%      Build Polynomial u(x,y)

% Build Polynomial u in x terms and y terms
u_x_poly = BuildPoly(u_x_roots);
u_y_poly = BuildPoly(u_y_roots);
% Build the polynomial u(x,y)
uxy_matrix_bi = u_x_poly' * u_y_poly;

% Strip the binomial coefficients from u_bi in the scaled Bernstein Basis,
% to obtain u in the standard Bernstein Basis
uxy_matrix_exact = zeros(m1-d1+1,m2-d2+1);
for i = 0:1:m1-d1
    uxy_matrix_exact(i+1,:) = uxy_matrix_bi(i+1,:) ./ nchoosek(m1-d1,i);
end
for j = 0:1:m2-d2
    uxy_matrix_exact(:,j+1) = uxy_matrix_exact(:,j+1) ./ nchoosek(m2-d2,j);
end


%%      Build Polynomial v(x,y)

% Build Polynomial v in x terms and y terms
v_x_poly = BuildPoly(v_x_roots);
v_y_poly = BuildPoly(v_y_roots);

% Build the polynomial v(x,y)
vxy_matrix_bi = v_x_poly' * v_y_poly;

% Strip the binomial coefficients from v_bi in the scaled Bernstein Basis,
% to obtain v in the standard Bernstein Basis
vxy_matrix_exact = zeros(n1-d1+1,n2-d2+1);
for i = 0:1:n1-d1
    vxy_matrix_exact(i+1,:) = vxy_matrix_bi(i+1,:) ./ nchoosek(n1-d1,i);
end
for j = 0:1:n2-d2
    vxy_matrix_exact(:,j+1) = vxy_matrix_exact(:,j+1) ./ nchoosek(n2-d2,j);
end


end





function u_x_roots = getQuotient(f_x_roots,d_x_roots)
%% given the roots of f in terms of x, and the roots of d in terms of x, 
% return the roots of the quotient polynomial u(x)

% Get number of roots in f(x)
num_roots_f_x = size(f_x_roots,1);

% Initialise an empty set of roots of u(x)
u_x_roots = [];

% for each of the roots in f(x)
for i = 1:1:num_roots_f_x
    % get the root
    root = f_x_roots(i,1);
    % get multiplicity
    mult_f = f_x_roots(i,2);
    
    
    % Look up in d_x
    if ~isempty(find(d_x_roots(:,1) == root))
        
        [row_d,~] = find(d_x_roots(:,1) == root);
        % if root is found, get multiplicity
        mult_d = d_x_roots(row_d,2);
        
        % subtract to obtain multiplicity in u(x)
        mult_u = mult_f - mult_d;
        if mult_u > 0
            u_x_roots = [u_x_roots; root mult_u];
        end
    else
        u_x_roots = [u_x_roots; root mult_f];
    end
end
end