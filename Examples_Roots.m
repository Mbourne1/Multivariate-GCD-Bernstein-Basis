function [fxy_matrix_exact,m] = Examples_Roots(ex_num)

% Each root is expressed in the form

%            ___
%   -r      |___|
%   1-r     |___|

    roots_f_x = [];
    roots_f_y = [];
    roots_f_xy = [];
    
switch ex_num
    
   
        
    
    case '1'
        % Relates to roots_examples tex file example 1
        
        % Get roots of f with respect to x
        roots_f_x =...
            [
            0.9   2;
            ]

        % Get roots of f with respect to y
        roots_f_y = ...
            [
            0.3   1
            ]

        % Get roots of f with respect to x and y
        % (x+y-0.1)
        roots_f_xy{1,1} = ...
            [
            0.1  1.1;
            1.1  2.1 
            ];
        
        m = 5;
        
    
        
    case '2'
        
        % Get roots and multiplicities wrt x
        roots_f_x = [...
            0.8     1;
            0.1     2;
            0.5     2;
            0.9     3;
            1.1234  4;
        ];
           
        m = 12;
        
    case '4'
        % Roots with respect to x
        roots_f_x = ...
            [
            0.5 1
            0.1 2
            ];
        % Roots with respect to y
        roots_f_y = ...
            [
            0.4 2
            ];
       
        m = 5;

 
    case '5'
        % (x-0.5)(x-0.1)^2 (x-0.7972)^3 (y-0.4)^2 (y-0.687) (x+y-0.1)
        
        % Roots with respect to x
        roots_f_x = ...
            [
            0.5     1
            0.1     2
            0.7972  3
            ];
        % Roots with respect to x
        roots_f_y = ...
            [
            0.4     2
            0.687   1
            ];
        % Roots with respect to x and y
        roots_f_xy{1,1} = ...
            [
            0.1  1.1;
            1.1  2.1 
            ];
        
        roots_f_xy{2,1} = ...
            [
            0.1  1.1;
            1.1  2.1 
            ];
        
     
        m = 13;
        case '5b'
        % (x-0.5)(x-0.1)^2 (x-0.7972)^3 (y-0.4)^2 (y-0.687) (x+y-0.1)
        
        % Roots with respect to x
        roots_f_x = ...
            [
            0.5     1
            0.1     2
            0.7972  3
            ];
        % Roots with respect to y
        roots_f_y = ...
            [
            0.4     2
            0.687   1
            ];
        % Roots with respect to x and y
        roots_f_xy{1,1} = ...
            [
            0.1  1.1;
            1.1  2.1 
            ];
        
        roots_f_xy{2,1} = ...
            [
            0.1  1.1;
            1.1  2.1 
            ];
         roots_f_xy{3,1} = ...
            [
            0.1  1.1;
            1.1  2.1 
            ];

        m = 15;
        
        case '6'
        
        % Roots with respect to x
        roots_f_x = ...
            [
            0.5     1
            0.1     2
            0.7972  3
            0.2475  4
            ];
        % Roots with respect to y
        roots_f_y = ...
            [
            0.4     2
            0.687   1   
            ];
        
        % Roots with respect to x and y
        roots_f_xy{1,1} = ...
            [
            0.1  1.1;
            1.1  2.1 
            ];
        
        roots_f_xy{2,1} = ...
            [
            0.1  1.1;
            1.1  2.1 
            ];
        roots_f_xy{3,1} = ...
            [
            0.1  1.1;
            1.1  2.1 
            ];
        
        m = 16;

        
end


roots_f_x = mult_roots_x(roots_f_x);
roots_f_y = mult_roots_y(roots_f_y);

roots_f = [roots_f_x; roots_f_y; roots_f_xy];

fxy_matrix_exact = BuildPoly_NonSeparable(roots_f);

end





