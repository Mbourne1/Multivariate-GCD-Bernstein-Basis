function [] = o_bivariate_mult(ex_num)


[f_x_roots,f_y_roots,...
    g_x_roots,g_y_roots,...
    u_x_roots, u_y_roots,...
    v_x_roots, v_y_roots,...
    d_x_roots, d_y_roots,...
    m1,m2,n1,n2,t1,t2] = examples(ex_num)

% 2.1 get coefficients of polynomial f

% % 2.1.1 get coefficients of polynomial f in terms of x
f_x_poly = BuildPoly(f_x_roots);

% % 2.1.2 get coefficients of polynomial f in terms of y
f_y_poly = BuildPoly(f_y_roots);

% % 2.1.3 Build matrix of coefficients of f(xy)
fxy_matrix = f_y_poly' * f_x_poly

fxy_vec = [];
for i = -m2 : 1 : m1
    i
    flipud(diag(flipud(fxy_matrix),i))
    fxy_vec = [fxy_vec ; flipud(diag(flipud(fxy_matrix),i))]
end

% 2.2 get coefficients of polynomial g

% % 2.2.1 get coefficients of polynomial g in terms of x
g_x_poly = BuildPoly(g_x_roots);

% % 2.2.2 get coefficients of polynomial g in terms of y
g_y_poly = BuildPoly(g_y_roots);

% % 2.2.3 get coefficients of polynomial g in terms of x and y
gxy_matrix = g_y_poly' * g_x_poly

gxy_vec = [];
for i = -n2 : 1 : n1
    gxy_vec = [gxy_vec ; flipud(diag(flipud(gxy_matrix),i))];
end

% 2.3 get coefficients of polynomial u

% % 2.3.1 get coefficients of polynomial u in terms of x
u_x_poly = BuildPoly(u_x_roots);

% % 2.3.2 get coefficients of polynomial u in terms of y
u_y_poly = BuildPoly(u_y_roots);

% % 2.3.3 get matrix of coefficients of u in terms of x and y
uxy_matrix = u_y_poly' * u_x_poly;

uxy_vec = []

try
    
    for i = -(m2-t2) : 1 : (m1-t1)
        i;
        flipud(diag(flipud(uxy_matrix),i))
        uxy_vec = [uxy_vec ; flipud(diag(flipud(uxy_matrix),i))];
    end
catch
    uxy_vec = uxy_matrix
end

uxy_vec

% 2.4 get coefficients of polynomial v

% % 2.4.1 get coefficients of polynomial v in terms of x
v_x_poly = BuildPoly(v_x_roots);

% % 2.4.2 get coefficients of polynomial v in terms of y
v_y_poly = BuildPoly(v_y_roots);

% % 2.4.3 get coefficients of polynomial v in terms of x and y
vxy_matrix = v_y_poly' * v_x_poly;

vxy_vec = [];


try
    for i = -(n2-t2) : 1 : (n1-t1)
        i;
        flipud(diag(flipud(vxy_matrix),i));
        vxy_vec = [vxy_vec ; flipud(diag(flipud(vxy_matrix),i))];
    end
catch
    vxy_vec = vxy_matrix;
end


% 2.5 get coefficients of polynomial f
% % 2.5.1 get coefficients of polynomial f in terms of x
d_x_poly = BuildPoly(d_x_roots);

% % 2.5.2 get coefficients of polynomial f in terms of y
d_y_poly = BuildPoly(d_y_roots);

% % 2.5.3 Build matrix of coefficients of f(xy)
dxy_matrix = d_y_poly' * d_x_poly

dxy_vec = [];
try
for i = -t2 : 1 : t1
    i
    flipud(diag(flipud(dxy_matrix),i))
    dxy_vec = [dxy_vec ; flipud(diag(flipud(dxy_matrix),i))]
end
catch
    dxy_vec = dxy_matrix';
end

% 3 perform multiplication between polynomials f and v

fxy_vec;
vxy_vec;

% 3.1 Build the Cauchy matrix for f

% get degree of u
% x component m1-t1
m1_t1 = m1-t1;
m2_t2 = m2-t2;

n1_t1 = n1-t1;
n2_t2 = n2-t2;

vxy_vec

% Build the fxy matrix surrounded by zeros such that the number of rows =
% m2+(n2-t)+1, and number of cols =
fxy_matrix


Cauchy_f = BuildCauchy(fxy_matrix,m1,m2,n1_t1,n2_t2);
Cauchy_g = BuildCauchy(gxy_matrix,n1,n2,m1_t1,m2_t2);

Cauchy_u = BuildCauchy(uxy_matrix,m1-t1,m2-t2,t1,t2);
Cauchy_v = BuildCauchy(vxy_matrix,n2-t1,n2-t2,t1,t2);

Cauchy_u * dxy_vec
fxy_vec
Cauchy_v * dxy_vec
gxy_vec

Cauchy_f * vxy_vec
Cauchy_g * uxy_vec


(Cauchy_f * vxy_vec)./(Cauchy_g * uxy_vec)




% number of columns in Cauchy matrix of f is given by m1-t1+1 x m2-t2+1


% 4 perform multiplication between polynomials g and u

% 4.1 Build the Cauchy matrix for g




end

function [f_bi] = BuildPoly(A)

% Calculate the number of distinct roots of the polynomial.
r = size(A,1);

% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.

f_bi = 1;
for k = 1:1:r
    w = B_conv(A(k,1),A(k,2));
    f_bi = conv(f_bi,w) ;
end


end

function [Cauchy_f] = BuildCauchy(fxy_matrix,m1,m2,n1_t1,n2_t2)


fxy_matrix_zeros = zeros(m2+(n2_t2)+1,m1+(n1_t1)+1);
fxy_matrix_zeros(1:size(fxy_matrix,1),1:size(fxy_matrix,2)) = fxy_matrix;

mycellarray = cell(n2_t2+1,n1_t1+1);
myvecarray = cell(n2_t2+1,n1_t1+1);

% for each column j in vxy matrix
for j = 0:1:n1_t1
    % for each row i in vxy matrix
    for i = 0:1:n2_t2
        temp_matrix = circshift(fxy_matrix_zeros,[i,j]);
        mycellarray{i+1,j+1} = temp_matrix;
        
        % Get vector of antidiagonals of temp matrix
        temp_vec = [];
        for k = -(m1+(n1_t1)+1) : 1 : (m2+(n2_t2)+1)
            temp_diag = flipud(diag(flipud(temp_matrix),k));
            temp_vec = [temp_vec ; temp_diag];
        end
        temp_vec;
        myvecarray{i+1,j+1} = temp_vec;
    end
end



Cauchy_f = [];

% for each matrix
% for each column
for k = 1:1:(n1_t1)+1
    % set the inital row to be one
    i=1;
    % decrease the index of the column by one
    for j = k:-1:1
        if i>(n2_t2)+1
            break;
        end
        myvecarray{i,j};
        Cauchy_f = [Cauchy_f myvecarray{i,j}];
        % increase the row index
        i = i+1;
    end
end

% for each row
% for i = 1...n2-t2+1, set j = max possible value then subtract one
for k = 2:1:n2_t2+1
    % Set the intial column to be the last column, and work to the first
    % column by subtracting one.
    j = n1_t1+1
    for i = k:1:n2_t2+1
        [i,j];
        Cauchy_f = [Cauchy_f myvecarray{i,j}];
        %decrease the column index
        j = j-1;
    end
end
end


function [t]=B_conv(rooot,mult)
%% This function convolves the vector [-r 1-r] with itself m times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% r : root
% m : multiplicity of root
% Outputs:
% t : vector which stores the result from this convolution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Note that (y-r) = -r(1-y) + (1-r)y and thus the polynomial y-r in the
% power basis is represented as the polynomial -r(1-y) + (1-r)y in the
% scaled Bernstein basis.


if mult==1
    t=[-rooot,1-rooot];
else
    
    q=[-rooot,1-rooot];
    t=[-rooot,1-rooot];
    for k=2:1:mult
        t=conv(t,q);
    end
end
end
