function [] = o_Bezier_Bezier_Intersection(ex_num)
% Given an example number, get the control points of 2 Bezier Curves, and
% calculate their intersections.
%
% 1. Get Control Points of polynomial f.
% 2. Get Control Points of polynomial g.
% 3. Implicitize Polynomial g.
% 4. Substitute the parametric expressions of f for x and y in t, into the
%    implicit polynomial g.

switch ex_num
    case '1'
        ex_num_f = '3';
        ex_num_g = '4';
end

% Get a set of control points for a polynomial f
CP_f = Examples_Bezier_ControlPoints(ex_num_f);

% Get the equation of the second bezier curve
CP_g = Examples_Bezier_ControlPoints(ex_num_g);


fprintf('Parametric representation of f in Brn Basis: \n')
CP_f;

% Get the vector of coefficients of f(x) and f(y)
f_x = CP_f(1,:);
f_y = CP_f(2,:);

f_x_Brn = Brn_Sym(CP_f(1,:));
f_y_Brn = Brn_Sym(CP_f(2,:));

fprintf('Parametric representation of g in Brn Basis: \n')
CP_g;
g_x_Brn = Brn_Sym(CP_g(1,:));
g_y_Brn = Brn_Sym(CP_g(2,:));

% Get the parametric expressions of polynomial f in the pwr basis
fprintf('Parametric representation of f in Power Basis:\n')
f_x_pwr_sym = Brn2Pwr_Sym(CP_f(1,:));
f_y_pwr_sym = Brn2Pwr_Sym(CP_f(2,:));

f_x_pwr = Bern2Power_Univariate(CP_f(1,:));
f_y_pwr = Bern2Power_Univariate(CP_f(2,:));

% Get the parametric expressions of polynomial g in the pwr basis
fprintf('Parametric representation of g in Power Basis:\n')
g_x_pwr = Brn2Pwr_Sym(CP_g(1,:));
g_y_pwr = Brn2Pwr_Sym(CP_g(2,:));

%%
% Implicitize the Bezier Control Points
fprintf('Implicit representation of g:\n')
[implicit_g,symbolic_expression] = Implicitize_Bezier_Sylvester(CP_g)

% Plot the implicit version of g against the parametric version of g

%%
% Get the parametric expressions of polynomial f in the brn basis
f_x =  CP_f(1,:)';
f_y =  CP_f(2,:);


% Put the sets of control points into an array
CP_arr{1} = CP_f;
CP_arr{2} = CP_g;

% Given the set of control points, draw the bezier curve
Plot_Bezier(CP_arr);

% Given the implicit representation of curve g. Plot it.
figure()
hold on
ezplot(symbolic_expression)
hold off



%% Substitute x(t) and y(t) into g

% get degree of input polynomial 
[~,c] = size(CP_f)
n = c-1;

[r,c] = size(implicit_g);

sum = zeros( (n^2)+1,1);

for i = 0:1:r-1
    for j = 0:1:c-1
        
               
        % Multiply x(t) by itself i times
        x_comp = 1;
        for k = 1:1:i   
            x_comp = Bernstein_Multiply(x_comp,f_x);
        end
        
        % Multiply y(t) by itself j times
        y_comp = 1;
        for k = 1:1:j
            y_comp = Bernstein_Multiply(y_comp,f_y');
        end
        
        %y_comp = y_comp';
        
        % Muliply x(t)^i and y(t)^j
        xy_comp = Bernstein_Multiply(x_comp,y_comp);
        
        uij =  xy_comp;
        
        % output polynomial will be of degree 2n
        m = n^2;
        % degree elevate uij
        [r1,~] = size(uij);
        curr_deg_uij = r1-1;
        num_deg_elv_req =  m-curr_deg_uij;

        uij = DegreeElevate_Univariate(uij,num_deg_elv_req);
        coef = implicit_g(i+1,j+1);
        uij = coef .* uij;
        
        sum = sum + uij
    end
end


end




% given the implicit matrix representation of


