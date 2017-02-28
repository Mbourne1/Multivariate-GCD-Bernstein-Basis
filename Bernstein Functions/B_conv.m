function [t] = B_conv(root,mult)
% B_conv(rooot,mult)
%
% This function convolves the vector [-r 1-r] with itself m times, where r
% is the root.
%
% Inputs:
%
%
% root : root
%
% mult : multiplicity of root
%
% Outputs:
%
%
% t : Vector which stores the coefficients of the polynomial obtained from
% convolution of the factor m times.



% Note that (y-r) = -r(1-y) + (1-r)y and thus the polynomial y-r in the
% power basis is represented as the polynomial -r(1-y) + (1-r)y in the
% scaled Bernstein basis.


if mult == 1 
    t=[-root;1-root];
else
    
    q =...
        [...
        -root;
        1-root
        ];
    t =...
        [...
        -root;
        1-root
        ];
    for k=2:1:mult
        t=conv(t,q);
    end
end