function arr_fww = GetPolynomialArrayWithThetas(arr_fxy, th1, th2)
% Given an array of polynomials, replace independent variables x_{1} =
% \theta_{1}\omega_{1} and x_{2} = \theta_{2}\omega_{2}


% Get number of polynomials

nPolys_arr_fxy = length(arr_fxy);

% Initialise array to store f_{i}(\omega_{1}, \omega_{2})
arr_fww = cell(nPolys_arr_fxy, 1);


for i = 1: 1 : nPolys_arr_fxy 
   
    % Get new polynomial f_{i}(\omega_{1}, \omega_{2})
    arr_fww{i} = GetWithThetas(arr_fxy{i}, th1, th2);
    
end


end