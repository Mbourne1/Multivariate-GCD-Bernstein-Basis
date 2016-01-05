function [hxy_matrix] = Bern_deconvolve_bivariate(fxy_matrix,gxy_matrix)

%%

% Get degrees of polynomial f
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get the degrees of polynomial g
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

% include the binomial coefficients in g
% Pre multiply gxy by the diagonal matrix nchoosek(n1,i)
Bi_n1 = zeros(1,n1+1);
for i = 0:1:n1
    Bi_n1(i+1) = nchoosek(n1,i);
end
Bi_n2 = zeros(1,n2+1);
for i = 0:1:n2
    Bi_n2(i+1) = nchoosek(n2,i);
end

% get g(x,y) including binomials
gxy_matrix_bi = diag(Bi_n1) * gxy_matrix * diag(Bi_n2);


% multiply gxy by the (m1-n1+1)*(m2-n2+1)

% build the matrix which will contain the multiplied entries
zero_matrix = zeros(m1+1,m2+1);

% the number of diagonals of the multiplication matrix h
num_diags = (m1 - n1 + 1) + (m2 - n2 + 1) + 1;

C1 = [];
% for every diagonal of the matrix dxy_mtrx.
for tot = 0:1:num_diags
    
    for i = tot:-1:0
        j = tot-i;

        % if i is within the number of rows of hxy_mtrx
        % if j is within the number of cols of hxy_mtrx
        if (i <= m1 - n1)  && (j <= m2 - n2)
            
            gxy_matrix_padded = zero_matrix;
            gxy_matrix_padded((i+1):(n1)+(i+1), (j+1):(n2)+(j+1)) = ...
                gxy_matrix_bi;
            
                        
            temp_vec = getAsVector(gxy_matrix_padded);
            
            C1 = [C1 temp_vec];
        end
    end
    
end

%% Get the polynomial f in vector form
f = getAsVector(fxy_matrix);

D = BuildD(0,0,m1,m2,0,0);
G = BuildG(m1-n1,m2-n2);
h = pinv(D*C1*G) * f;


% Perform QR decomposition to obtain appoximation x_ls 
DCQ = D*C1*G;
bk = f;
[~,n3] = size(DCQ);
[Q1,R] = qr(DCQ);

R1 = R(1:n3,:);
cd = Q1'*bk;
c = cd(1:n3,:);
x_ls = R1\c;

% set h1 to the least squares solution.
h1 = x_ls;

hxy_matrix = getAsMatrix(h,m1-n1,m2-n2);
h1_matrix = getAsMatrix(h1,m1-n1,m2-n2);
hxy_matrix;
end