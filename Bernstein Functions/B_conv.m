function [t] = B_conv(root, multiplicity)
% B_conv(rooot, mult)
% 
% This function convolves the vector [-r 1-r] with itself m times, where r
% is the root.
%
% Inputs:
%
%
% root : (Float) Root
%
% multiplicity : (Int) Multiplicity of root
%
% Outputs:
%
%
% t : (Vector) Coefficients of the polynomial obtained from convolution of 
% the factor (x - r) m times.



% Note that (y-r) = -r(1-y) + (1-r)y and thus the polynomial y-r in the
% power basis is represented as the polynomial -r(1-y) + (1-r)y in the
% scaled Bernstein basis.


if multiplicity == 1 
    t = [-root; 1 - root];
else
    
    q =...
        [...
        -root;
        1 - root
        ];
    t =...
        [...
        -root;
        1 - root
        ];
    
    for k = 2 : 1 : multiplicity
    
        t=conv(t,q);
    end
    
end