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
partial_fxy = zeros(m1+1,m2);

% for each row of coefficients
for i = 0:1:m1
    curr_row = fxy(i+1,:);
    new_row = diff_y(curr_row);
    partial_fxy(i+1,:) = new_row;
end


end

function fxy_wrt_y = diff_y(fxy)
% Differentiate the given vector with respect to x.


% get the size of polynomial fxy
[~,m2] = GetDegree_Bivariate(fxy);


mult_mat = zeros(m2+1,m2);

for i = 0:1:m2-1
    mult_mat(i+1,i+1) = -m2;
    mult_mat(i+2,i+1) = m2;
end

fxy_wrt_y = fxy * mult_mat  ;

fxy_wrt_y = fxy_wrt_y;
end