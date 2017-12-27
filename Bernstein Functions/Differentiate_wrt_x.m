function  [partial_fxy] = Differentiate_wrt_x(fxy)
% Given the matrix of coefficients of the bivariate polynomial f(x,y),
% diffferentiate with respect to x.
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% % Outputs
%
% partial_fxy : (Matrix) Coefficients of the partial derivative of f(x,y)


% Get degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Initialise a matrix to store the partial derivative of f(x,y) with
% respect to x
partial_fxy = zeros(m1, m2 + 1);

% for each column of coefficients
for j = 0:1:m2
    
    % Get Column from f(x,y)
    current_column = fxy(:, j + 1);
    
    % Differentiate the column with respect to x
    new_column = diff_x(current_column);
    
    % Place differentiated column into partial f(x,y) coefficient matrix
    partial_fxy(:, j + 1) = new_column;
    
end


end

function fxy_wrt_x = diff_x(fxy)
% Differentiate the given vector with respect to x.
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% % Outputs
% 
% fxy_wrt_x : (Matrix) Coefficients of partial derivative of f(x,y) with 
% respect to x 



% Get degree of f(x,y) with respect to x.
[m1, ~] = GetDegree_Bivariate(fxy);

% Transpose f(x,y)
fxy = fxy';

mult_mat = zeros(m1 + 1, m1);

for i = 0 : 1 : m1 - 1
    mult_mat(i + 1, i + 1) = -m1;
    mult_mat(i + 2, i + 1) = m1;
end

fxy_wrt_x = fxy * mult_mat;

fxy_wrt_x = fxy_wrt_x';
end