function cellArr = mult_roots_y(root_mult_mat)
%
% Given the root and multiplicity matrix produce a cell array of the roots
% where multiple roots r_{i} are repeated m_{i} times, where m_{i} is the 
% multiplicity of r_{i}.
%
% % Inputs:
%
%                    _______________
% root_mult_mat =   | root  |  mult |
%                   |       |       |
%                   |_______|_______|
%
% % Outputs
%
% cellArr : (Array) 

%
% Initialise a count
count = 1;

% Initialise a cell array
cellArr = {};

% Get number of distinct roots
[nDistinctFactors] = size(root_mult_mat,1);

% For each distinct root in the array
for i = 1:1:nDistinctFactors
    
    % Get the ith root
    root = root_mult_mat(i,1);
    
    % Get the multiplicity of the ith root
    m = root_mult_mat(i,2);
    
       
    % Insert root into cell array m times
    for j = 1:1:m
        cellArr{count, 1} = root_y(root);
        count = count + 1;
    end
end
end