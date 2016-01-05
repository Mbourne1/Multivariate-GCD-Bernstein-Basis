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
        % Get roots with respect to x
        roots_f_x = ...
            [
            0.5 1
            0.1 2
            ]
        roots_f_y = ...
            [
            0.4 2
            ]
       
        m = 5;

 
    case '5'
        % (x-0.5)(x-0.1)^2 (x-0.7972)^3 (y-0.4)^2 (y-0.687) (x+y-0.1)
        
        % Get roots with respect to x
        roots_f_x = ...
            [
            0.5     1
            0.1     2
            0.7972  3
            ];
        roots_f_y = ...
            [
            0.4     2
            0.687   1
            ];
        
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
        
        % Get roots with respect to x
        roots_f_x = ...
            [
            0.5     1
            0.1     2
            0.7972  3
            ];
        roots_f_y = ...
            [
            0.4     2
            0.687   1
            ];
        
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
        
        % Get roots with respect to x
        roots_f_x = ...
            [
            0.5     1
            0.1     2
            0.7972  3
            0.2475  4
            ]
        roots_f_y = ...
            [
            0.4     2
            0.687   1   
            ]
        
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
        
        roots_f_x = mult_roots_x(roots_f_x);
        roots_f_y = mult_roots_y(roots_f_y);
        
        roots_f = [roots_f_x ; roots_f_y ; roots_f_xy];
        fxy_matrix_exact = BuildPoly_NonSeparable(roots_f);
        m = 16;
%     case '1'
%         
%         % Root of f in terms of x
%         roots_f{1,1} = ...
%             [
%             -0.7;
%             0.3
%             ];
%         
%         fxy_matrix_exact = BuildPoly_NonSeparable(roots_f)
%         
%         m = 1;
%         
%     case '2'
%         
%         % Roots of f in terms of x
%         roots_f{1,1} = ...
%             [
%             -0.8
%             0.2
%             ];
%         roots_f{2,1} = ...
%             [
%             -0.8
%             0.2
%             ];
%         roots_f{3,1} = ...
%             [
%             -0.3
%             0.7
%             ];
%         roots_f{4,1} = ...
%             [
%             -0.2     0.8
%             ];
%         roots_f{5,1} = ...
%             [
%             -0.2     0.8
%             ];
%         roots_f{6,1} = ...
%             [
%             0.1    1.1     
%             ];
%         
%         fxy_matrix_exact = BuildPoly_NonSeparable(roots_f)
%         m = 6;
%         
%     case '3'
%         roots_f{1,1} = ...
%             [
%             -1  0
%             0   1
%             ];
%         roots_f{2,1} = ...
%             [
%             -1;
%             0
%             ];
%         roots_f{3,1} = ...
%             [
%             -5  -4
%             ];
%         fxy_matrix_exact = BuildPoly_NonSeparable(roots_f)
%         m = 3;
%     
%     case '30'
%         roots_f{1,1} = ...
%             [
%             -5  -2      0;
%             0   1.25    2;
%             0   0       0
%             ]
%         fxy_matrix_exact = BuildPoly_NonSeparable(roots_f)
%         m = 4;
%         
%     case '31'
%        roots_f{1,1} = ...
%            [
%             600 1280/3  1398/5 831/5    256/3   32  0
%             250 1450/9  1585/18 343/10  -62/45  -200/9  -32
%             70  298/9   902/225 -4753/300   -6134/225   -32 -32 
%             45/4    -3/4    -1463/150   -3049/200   -1321/75    -88/5 -16
%             50  95/2    686/15  2241/50 10046/225   2026/45 686/15
%             250 250 250 250 250 250 250
%             750 750 750 750 750 750 750
%             ]
%         fxy_matrix_exact = BuildPoly_NonSeparable(roots_f)
%         m = 12;
%         
%     case '32'
%        roots_f{1,1} = ...
%            [
%             -4  -4  -3
%             -4  -4  -3  
%             -3  -3  -2
%            ]
%        
%        
%         fxy_matrix_exact = BuildPoly_NonSeparable(roots_f)
%        m = 12;
        
end


roots_f_x = mult_roots_x(roots_f_x);
roots_f_y = mult_roots_y(roots_f_y);

roots_f = [roots_f_x; roots_f_y; roots_f_xy];

fxy_matrix_exact = BuildPoly_NonSeparable(roots_f);

end





