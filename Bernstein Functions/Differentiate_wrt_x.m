function  [partial_fxy] = Differentiate_wrt_x(fxy)
% Given the matrix of coefficients of the bivariate polynomial f(x,y),
% diffferentiate with respect to x.
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)

% Get degree of f(x,y)
[m1,m2] = GetDegree_Bivariate(fxy);

% Initialise a matrix to store the partial derivative of f(x,y) with
% respect to x
partial_fxy = zeros(m1,m2+1);

% for each column of coefficients
for j = 0:1:m2
    
    curr_col = fxy(:,j+1);
    
    new_col = diff_x(curr_col);
    
    partial_fxy(:,j+1) = new_col;
    
end


end

function fxy_wrt_x = diff_x(fxy)
% Differentiate the given vector with respect to x.
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)

% Get degree of f(x,y) with respect to x.
[m1,~] = GetDegree(fxy);

% Transpose f(x,y)
fxy = fxy';

mult_mat = zeros(m1+1,m1);

for i = 0:1:m1-1
    mult_mat(i+1,i+1) = -m1;
    mult_mat(i+2,i+1) = m1;
end

fxy_wrt_x = fxy * mult_mat;

fxy_wrt_x = fxy_wrt_x';
end