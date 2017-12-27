function  [partial_fxy] = Differentiate_wrt_y(fxy)
% Given the matrix of coefficients of the bivariate polynomial f(x,y),
% diffferentiate with respect to y
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% % Outputs
%
% partial_fxy : (Matrix) Coefficients of the derivative of f(x,y)

% Get the degree of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Initialise matrix
partial_fxy = zeros(m1 + 1, m2);

% for each row of coefficients of f(x,y)
for i = 0 : 1 : m1

    % Get the contents of the current row
    curr_row = fxy(i + 1,:);

    % Get derivative of current row with respect to y
    new_row = diff_y(curr_row);

    % Insert contents into matrix of coefficients of partial derivative
    partial_fxy(i+1,:) = new_row;

end


end

function fxy_wrt_y = diff_y(fxy_row)
% Differentiate the given vector with respect to x.
%
% % Inputs
%
% fxy_row : (Vector) Coefficients of one row of polynomial f(x,y)
%
% % Outputs
%
% fxy_wrt_y : (Vector) Coefficients of row of partial derivative of f(x,y)


% Get the degree of the polynomial with respect to y
[~, m2] = GetDegree_Bivariate(fxy_row);

% Create the transformation matrix
mult_mat = zeros(m2 + 1, m2);

for i = 0 : 1 : m2 - 1

    mult_mat(i + 1, i+1) = -m2;
    mult_mat(i + 2, i+1) = m2;

end

% Get parital derivative
fxy_wrt_y = fxy_row * mult_mat  ;


end