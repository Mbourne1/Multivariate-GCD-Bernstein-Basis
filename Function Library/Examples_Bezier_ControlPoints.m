function CP = Examples_Bezier_ControlPoints(ex_num)
switch ex_num
    
    case '1'
        CP = ...
        [
            0.1     0.5     0.8 
            0.5     2       1
        ];
    case '2'

        CP = ...
        [
            0.1     0.5     0.8
            1       -0.7    1.1
        ];
    
    
    case '3'
        
        CP = ...
        [
            0.1     0.5     0.8     1.2
            1       -0.7    1.1     0.5
        ];
    case '4'
        CP = ...
            [
             0.1     0.5     0.8     1.2
            1       -0.8    1.5     0.5
            ];
        

end