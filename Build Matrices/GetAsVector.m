function v_fxy = GetAsVector(fxy)
% Given the polynomial f(x,y) in the Bernstein basis, whose coefficients are in
% matrix form, obtain the vector of the coefficients such that the order is
% increasing and the higher power of x is first.
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) as a matrix
%
% % Outputs
%
% f_vec : (Vector) Vector of ordered coefficients of f(x,y)

global SETTINGS


switch SETTINGS.VECTORISATION_METHOD
    
    case 'Version 1'

        v_fxy = GetAsVector_Version1(fxy);

    case 'Version 2'

        v_fxy = GetAsVector_Version2(fxy);

    otherwise

        error('Error')
end

end