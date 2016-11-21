function [] = o_GetIntersection_Parametric_Parametric(ex_num)
% Given two parametrically defined curves, get the intersection points.


num_data_points = 1000;
t= linspace (-1, 2, num_data_points);

switch ex_num
    
        
    case '2'
        % The parametrically defined curve
        
        
        X1 = (-0.2).*(t.^2) + (0.8.*t) + 0.1;
        Y1 = (-2.5).*(t.^2) + (  3.*t) + 0.5;
        
        % The implicitly defined curve
        X2 = t;
        Y2 = 2.*t;
        
        X3 = t;
        Y3 = -2.1.*(t.^2) + 1.4.*t +0.3
        
    case '3'
        X1 = 0.1 + 0.4.*t
        Y1 = 0.8 + 0.4.*t;
        
        X2 = t;
        Y2 = -t + 1.3;
    case '4'
         % Intersection of two bezier curves
         X1 = -0.1*(t.^2) + 0.8*t + 0.1;
         Y1 = -2.5*(t.^2) + 3*t + 0.5;
        
         X2 = -0.1.*(t.^2) + 0.8*t +0.1;
         Y2 =  3.5.*(t.^2) - 3.4*t + 1;
         
end

parametric{1} = [X1; Y1];
parametric{2} = [X2; Y2];


Plot_Parametric(parametric)



end