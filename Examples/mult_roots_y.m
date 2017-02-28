
function cellArr = mult_roots_y(root_mult_mat)
%%
% Given the root and multiplicity matrix produce a cell array of the roots
% where multiple roots r_{i} are repeated m_{i} times, where m_{i} is the 
% multiplicity of r_{i}.

% Inputs:

%                    _______________
% root_mult_mat =   | root  |  mult |
%                   |       |       |
%                   |_______|_______|

%%
% Initialise a count
count = 1;

% Initialise a cell array
cellArr = {};

% Get number of distinct roots
[num_roots,~ ] = size(root_mult_mat);

% For each distinct root in the array
for i = 1:1:num_roots
    
    % Get the multiplicity of the root
    m = root_mult_mat(i,2);
    
    % Get the root
    root = root_mult_mat(i,1);
    
    % Insert root into cell array m times
    for j = 1:1:m
        cellArr{count,1} = root_y(root);
        count = count + 1;
    end
end
end