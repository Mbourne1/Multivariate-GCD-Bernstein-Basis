function arr_hxy = Deconvolve_Bivariate_Separate(arr_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices)
%
% % Outputs
%
% arr_hxy : (Array of Matrices)



% Get number of polynonmials in array
nPolys_fxy = length(arr_fxy);
nPolys_hxy = nPolys_fxy - 1;

% Initialise an array to store h_{i}(x,y)
arr_hxy = cell(nPolys_hxy,1);


for i = 1 : 1 : nPolys_fxy - 1

    fxy = arr_fxy{i};
    gxy = arr_fxy{i + 1};
    
    arr_hxy{i} = Deconvolve_Bivariate(fxy,gxy);
    
end



end