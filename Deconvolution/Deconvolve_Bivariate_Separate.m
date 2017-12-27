function arr_hxy = Deconvolve_Bivariate_Separate(arr_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices) Each matrix contains coefficients of a
% polynomial f_{i}(x, y)
%
% % Outputs
%
% arr_hxy : (Array of Matrices) Each matrix contains coefficients of a
% polynomail h_{i}(x, y)



% Get number of polynomials in array f_{i}(x, y)
nPolys_fxy = length(arr_fxy);

% Get number of polynomials in array h_{i}(x, y)
nPolys_hxy = nPolys_fxy - 1;

% Initialise an array to store h_{i}(x,y)
arr_hxy = cell(nPolys_hxy, 1);

% Perform set of (n - 1) deconvolutions
for i = 1 : 1 : nPolys_fxy - 1

    
    fxy = arr_fxy{i};
    gxy = arr_fxy{i + 1};
    
    % Store result of deconvolution in array h_{i}(x,y)
    arr_hxy{i} = Deconvolve_Bivariate(fxy,gxy);
    
end



end