function [] = o_GetIntersection_Implic_Parametric()
% given an implicitly defined curve, and a parametrically defined curve.
% get the intersection points.

ex_num = '2';
num_data_points = 100;
t= linspace (-4, 4, num_data_points);

pwr_basis = [t.^2; t ; ones(1,num_data_points)];

switch ex_num
    case '1'
        % The parametrically defined curve
        X = ((3*t)./(1+t.^3));
        Y = ((3*t.^2)./(1+t.^3));
        
    case '2'
        % The parametrically defined curve in power basis
        
        % The coefficients of parametrically defined x(t) where leading
        % coefficients has highest power.
        coeff_X =  [-0.2   0.8   0.1]; 
        coeff_Y =  [-2.5   3     0.5]; 
        
        % Get the vector pwr_basis so that when multiplied we generate the
        % series of datapoints
        degree_f = size(coeff_X,2)        
       
        X = coeff_X  * pwr_basis;
        Y = coeff_Y  * pwr_basis;
        
        % The implicitly defined curve
        X2 = t;
        Y2 = 2*t;
        
        % Substitute the parametric form into the implicit form
        X3 = t;
        Y3 = 2.1*(t.^2) - 1.4*t -0.3;
        
end

parametric = [X; Y];
implicit = [X2; Y2];






end