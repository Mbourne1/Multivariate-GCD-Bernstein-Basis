function arr_fxy = GetPolynomialArrayFromVector(v_fxy, vDeg_x_fxy, vDeg_y_fxy)
% % Inputs
% 
% v_fxy : (Vector) Vector containing coefficients of all polynomials f_{i}(x,y)
%
% vDeg_x_fxy : (Vector) Vector containing the degree of each f_{i}(x,y)
% with respect to x
%
%
% vDeg_y_fxy : (Vector) Vector containing the degree of each f_{i}(x,y)
% with respect to y
%
% % Outputs
%
% arr_fxy : (Array of Matrices)


nPolys_arr_fxy = size(vDeg_x_fxy,1);
arr_fxy = cell(nPolys_arr_fxy, 1);


for i = 1 : 1 : nPolys_arr_fxy
   
    % Get degree of f_{i}(x,y)
    m1 = vDeg_x_fxy(i);
    m2 = vDeg_y_fxy(i);
    
    % Get number of coefficients
    nCoefficients_fxy = (m1 + 1) * (m2 + 1);
    
    % Get temporary vector containing only coefficients of f_{i}(x,y)
    v_temp = v_fxy(1 : nCoefficients_fxy);
    
    % Get matrix of coefficients 
    arr_fxy{i} = GetAsMatrix(v_temp, m1, m2);
    
    % Remove coefficients from vector
    v_fxy(1:nCoefficients_fxy) = [];
    
end


end