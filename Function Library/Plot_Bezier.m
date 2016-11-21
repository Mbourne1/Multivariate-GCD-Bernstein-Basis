function [] = Plot_Bezier(cp_f_arr)
% Given a set of control points for curves {C1, C2,...} plot the curves
% Input is an array of polynomial control points


% Plot the Graph
figure('name','Bezier Plot')
hold on

% Get the number of sets of control points
[r,c] = size(cp_f_arr);

for curve_num = 1:1:c
    
    % Get the current set of control points
    cp_f = cp_f_arr{curve_num};
    
    % Given a set of control points, plot the bezier curve.
    degree_f = size(cp_f,2) - 1;
    
    t = linspace(0,1,101);
    
    pts_f = 0;
    for i = 0:1:degree_f
        
        pts_f = pts_f + kron( nchoosek(degree_f,i).*((1-t).^(degree_f - i)) .* (t.^(i)) ,cp_f(:,i+1));
    end

    % For every column (control point) of f
    for i = 1:1:degree_f+1
        str = sprintf('c_%i',i);
        placelabel(cp_f(:,i),str);
    end
    
    plot_name = sprintf('Curve C_%i',curve_num);
    plot(pts_f(1,:),pts_f(2,:),'DisplayName',plot_name)
    
     
       
end
legend('show')
hold off

end

function placelabel(pt,str)
x = pt(1);
y = pt(2);
h = line(x,y);
h.Marker = '.';
h = text(x,y,str);
h.HorizontalAlignment = 'center';
h.VerticalAlignment = 'bottom';
end