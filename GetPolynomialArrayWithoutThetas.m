function arr_fxy = GetPolynomialArrayWithoutThetas(arr_fww, th1, th2)
%
% % Inputs
%
% arr_fww : (Array of Matrices)
%
% th1 : (Float)
%
% th2 : (Float)
%
% % Outputs
%
% arr_fxy : (Array of Matrices)


% Get number of polynomials in the array 
nPolys_fxy = size(arr_fww,1);

% Initialise array to store f_{i}(x,y)
arr_fxy = cell(nPolys_fxy,1);

for i = 1 : 1 : nPolys_fxy
    
    % Get f_{i}(x,y)   
    arr_fxy{i} = GetWithoutThetas(arr_fww{i}, th1, th2);
    
end




end