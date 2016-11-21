function [A] = Examples_Roots_fx(ex_num)
% Given an example number, return the matrix whose rows row_{i} consist of
% a root r_{i} and its multiplicity m_{i}.
%
% % Input.
%
% ex_num :
%
% % Output:
%
% A : Matrix where the first column contains the roots r_{i}, and the
% second column contains the multiplicities m_{i} of the roots.

switch ex_num
    
    case '1'
        
        A = [
            2   1
            3   2
            ];
    case '2'
        A = [
            1   1
            2   2
            3   1
            ];
    case '3'
        A = [
            1 1
            2 2
            3 3
            ];
    otherwise
        error('Example corresponding to the input example number does not exist.')
end


end
