function [] = PrintRoots_Bivariate(root_mult_arr_x, root_mult_arr_y, polys_xy)
% PrintRoots_Bivariate(root_mult_arr_x, root_mult_arr_y, polys_xy)
%
% Given the roots of f(x,y) and their corresponding multiplicities, print
% out the factorisation of f(x,y).
%
% Inputs.
%
%
% root_mult_arr_x : A matrix of roots and multiplicities in terms of x.
%
% root_mult_arr_y : A matrix of roots and multiplicities in terms of y.
%
% polys_xy : Array of bivariate polynomials which are the non-separable
%            factors of f(x,y)

nRoots_x = size(root_mult_arr_x, 1);
nRoots_y = size(root_mult_arr_y, 1);

nFactors_xy = size(polys_xy,1);

% for each root in terms of x
str = sprintf('f(x) =' );

str = strcat(str,'... \n');

for i = 1:1:nRoots_x
    
    % Get root.
    root = root_mult_arr_x(i,1);
    
    % Get Multiplicity of root
    mult = root_mult_arr_x(i,2);
    
    % Get the first string
    str_a = '(x ';
    
    % Get the sign string
    if root >=0
        
        str_sign = '-';
        
    else % Root is negative
        
        str_sign = '+';
    end
    
    str_b = sprintf('%2.4f)^{%i}', abs(root), mult);
    
    temp_str = strcat(str_a,str_sign,str_b);
    
    str = strcat(str,temp_str);
end

str = strcat(str,'... \n');

for i = 1:1:nRoots_y
    
    root = root_mult_arr_y(i,1);
    mult = root_mult_arr_y(i,2);
    
    str_a = '(y ';
    
    if root >=0
        
        str_sign = '-';
        
    else % Root is negative
        
        str_sign = '+';
    end
    
    str_c = sprintf('%2.4f)^{%i}', abs(root), mult);
    
    temp_str = strcat(str_a,str_sign,str_c);
    
    str = strcat(str,temp_str);
end

str = strcat(str,'... \n');

for k = 1:1:nFactors_xy
    
    coef_matrix = polys_xy{k};
    
    [m1,m2] = GetDegree_Bivariate(coef_matrix);
    
    % Get the number of diagonals in the matrix of coefficients of f(x,y)
    nDiags_fx = (m1+1)+(m2+1)-1;
    
    temp_str2 = '(';
    
    
    for tot = 0:1:nDiags_fx
        
        
        
        for i = tot:-1:0
            j = tot - i;
            
            if(i > m1 || j > m2)
                % Do Nothing
            else
                
                % Get the coefficient
                coefficient = coef_matrix(i+1,j+1);
                
                % Get Sign string
                if coefficient > 0
                    if i == 0 && j == 0
                        % If first coeffficient, and sign is positive, then
                        % drop the sign
                        str_sign = '';
                    else
                        % If not the first coefficient
                        str_sign = '+';
                    end
                else
                    % Negative coefficient
                    str_sign = '-';
                end
                
                
                % Get the coefficient string
                str_coef = sprintf('%2.2f',abs(coefficient));
                
                % Get the x power string
                if i ~= 0
                    str_x = sprintf('x^{%i}',i);
                else
                    str_x = '';
                end
                
                % Get the y power string
                if j ~= 0
                    str_y = sprintf('y^{%i}',j);
                else
                    str_y = '';
                end
                
                % Create the string which has sign, coefficient and variables
                temp_str = strcat(str_sign,str_coef,str_x,str_y);
                
                % Add the string to the tempo string
                temp_str2 = strcat(temp_str2,temp_str);
            end
            
            
            
            
            
        end
        
    end
    
    % Close bracket
    temp_str2 = strcat(temp_str2, ')');
    
    str = strcat(str,temp_str2);
    
    fprintf(str)
    fprintf('\n')
    
    
    
    
    
    
end