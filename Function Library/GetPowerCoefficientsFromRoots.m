function [fx] = GetPowerCoefficientsFromRoots(ex_num)
% Given a set of polynomial roots and multiplicities, construct the power 
% basis coefficients of the polynomial.


A = Examples_Roots_fx(ex_num)


% Get the coefficients of polynomial f in terms of power basis, where the
% leading coefficient is of highest degree.
fx = P_poly(A);



end

function [f_bi] = P_poly(A)
% Obtain polynomial coefficients in the Power Basis, given a set of roots 
% and multiplicities. 
% This function implements the convolution operation on the polynomial
% defined by the matrix A, which has two columns, that is, this function
% computes the coefficients of the polynomial defined by A. These are
% the coefficients in the scaled Bernstein basis form of the polynomial.
% Coefficients of the form ai(m choose i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the number of distinct roots of the polynomial.
    r = size(A,1); 

% Convolve each factor, which is defined by a row of A, separately. 
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.
  
    f_bi = 1;
    for k = 1:1:r
        w = P_conv(A(k,1),A(k,2));    
        f_bi = conv(f_bi,w) ;
    end    
    f_bi = f_bi';
end


function [t]=P_conv(root,mult)
% This function convolves the vector [-r 1-r] with itself m times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% r : root
% m : multiplicity of root
% Outputs:
% t : vector which stores the result from this convolution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Note that (y-r) = -r(1-y) + (1-r)y and thus the polynomial y-r in the
% power basis is represented as the polynomial -r(1-y) + (1-r)y in the 
% scaled Bernstein basis.   


    if mult==1
        t=[1,-root];
    else

        q=[1,-root];
        t=[1,-root];
        for k=2:1:mult
            t=conv(t,q);
        end
    end
end