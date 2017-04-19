function fxy_vec = GetAsVector_Version2(fxy)
% Given the coefficient matrix of the polynomial f(x,y), get the
% coefficients as a vector
% 
% Coefficients are ordered based on columns of matrix. 
% f_vec = [c1 ; c2 ; ...; c_{m2}] where c_{i} are the m_{2}+1 columns of
% the coefficient matrix
%
% % Inputs
%
% fxy : (Matrix) Matrix of coefficients of f(x,y)
%
% % Outputs
%
% fxy_vec : (Vector) Ordered Coefficients of f(x,y)

% Get the degree structure of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get number of coefficients in the polynomial f(x,y)
nCoefficients = (m1+1) * (m2+1);

% Initialise the matrix which stores the coefficients
fxy_vec = zeros(nCoefficients, 1);


% Build the matrix of coefficients
for j = 0 : 1 : m2
    
    start_point = (j*(m1+1)) + 1;
    
    temp_vec = fxy(:,j+1);
    
    fxy_vec(start_point:start_point + (m1),1) = temp_vec;
    
end



end