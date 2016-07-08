function [fxy,gxy,dxy,m,n,t] = Examples_GCD(ex_num)

EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
    case 'From Roots'
        [fxy,gxy,dxy,m,n,t] = Examples_GCD_FromRoots(ex_num);
        
    case 'From Coefficients'
        [fxy,gxy,dxy,m,n,t] = Examples_GCD_FromCoefficients(ex_num);
        
    otherwise
        error('err');
end


end

