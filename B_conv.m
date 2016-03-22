function [t] = B_conv(rooot,mult)
% B_conv(rooot,mult)
% This function convolves the vector [-r 1-r] with itself m times, where r
% is the root.
%
% Inputs:
%
%
% r : root
%
% m : multiplicity of root
%
% Outputs:
%
%
% t : Vector which stores the result from this convolution.
%


% Note that (y-r) = -r(1-y) + (1-r)y and thus the polynomial y-r in the
% power basis is represented as the polynomial -r(1-y) + (1-r)y in the
% scaled Bernstein basis.


if mult==1
    t=[-rooot;1-rooot];
else
    
    q =...
        [...
        -rooot;
        1-rooot
        ];
    t =...
        [...
        -rooot;
        1-rooot
        ];
    for k=2:1:mult
        t=conv(t,q);
    end
end