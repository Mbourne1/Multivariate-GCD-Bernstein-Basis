
function [vRHS] = BuildRHS_vec(arr_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices)
%
% % Outputs
%
% vRHS


nPolys_arr_fxy = length(arr_fxy);

vRHS = [];

for i = 1 : 1 : nPolys_arr_fxy - 1
    
    fxy = arr_fxy{i};
    
    temp_vec = GetAsVector(fxy);
    
    vRHS = [vRHS ; temp_vec];
    
end





end