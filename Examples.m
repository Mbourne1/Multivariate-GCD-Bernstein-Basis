function [f_x_roots, f_y_roots, g_x_roots, g_y_roots, d_x_roots, ...
    d_y_roots,u_x_roots,u_y_roots, v_x_roots, v_y_roots] = Examples(ex_num)

switch ex_num
    
    case -1
        f_x_roots = [
            2   1;
            3   1;
            ];
        f_y_roots = [
            2   1;
            5   1;
            ];
        
        g_x_roots = [
            2   1;
            1   1;
            ];
        g_y_roots = [
            2   1;
            6   1;
            ];
        
        d_x_roots = [
            2   1;
            ];
        
        d_y_roots = [
            2   1;
            ];
        
        
    case 0
        f_x_roots = [
            0.2   1;
            0.3   1;
            ];
        f_y_roots = [
            0.2   1;
            0.3   1;
            ];
        
        g_x_roots = [
            0.2   1;
            0.7   1;
            ];
        g_y_roots = [
            0.2   1;
            0.7   1;
            ];
        
        d_x_roots = [
            0.2   1;
            ];
        
        d_y_roots = [
            0.2   1;
            ];
        
    case 1
        f_x_roots = [
            0.2   1;
            0.3   1;
            0.4   1;
            0.9   1;
            ];
        f_y_roots = [
            0.2   1;
            0.3   1;
            0.4   1;
            0.5   1;
            0.6   1;
            ];
        
        g_x_roots = [
            0.2   1;
            0.4   1;
            0.7   1;
            0.9   1;
            ];
        g_y_roots = [
            0.2   1;
            0.5   1;
            0.8   1;
            ];
        
        d_x_roots = [
            0.2   1;
            0.4   1;
            0.9   1;
            ];
        
        d_y_roots = [
            0.2   1;
            0.5   1;
            ];
        
    
      
        
    case 3
        f_x_roots = [
            0.3   1;
            0.4   1;
            0.9   1;
            ];
        f_y_roots = [
            0.2   1;
            0.3   1;
            0.5   1;
            ];
        
        g_x_roots = [
            0.3   1;
            0.4   1;
            ];
        
        g_y_roots = [
            0.1   1;
            0.5   1;
            0.8   1;
            ];
        
        d_x_roots = [
            0.3   1;
            0.4   1;
            ];
        
        d_y_roots = [
            0.5   1;
            ];
        
    case 4
        f_x_roots = [
            0.3   2;
            0.4   3;
            0.9   1;
            ];
        f_y_roots = [
            0.2   1;
            0.3   1;
            0.5   1;
            ];
        
        g_x_roots = [
            0.3   2;
            0.4   3;
            ];
        
        g_y_roots = [
            0.1   1;
            0.5   1;
            0.8   1;
            ];
        
        d_x_roots = [
            0.3   2;
            0.4   3;
            ];
        
        d_y_roots = [
            0.5   1;
            ];
        
    case 5
        f_x_roots = [
            0.3   2;
            0.4   3;
            0.9   1;
            0.1   1;
            ];
        f_y_roots = [
            0.2   1;
            0.3   1;
            0.5   1;
            ];
        
        g_x_roots = [
            0.3   2;
            0.4   3;
            0.1   1;
            ];
        
        g_y_roots = [
            0.1   1;
            0.5   1;
            0.8   1;
            ];
        
        d_x_roots = [
            0.3   2;
            0.4   3;
            0.1   1;
            ];
        
        d_y_roots = [
            0.5   1;
            ];
        
        case 6
        f_x_roots = [
            0.3     2;
            0.4     3;
            0.9     1;
            0.1     1;
            1.1     1;
            ];
        f_y_roots = [
            0.2   1;
            0.3   1;
            0.5   1;
            ];
        
        g_x_roots = [
            0.3   2;
            0.4   3;
            0.1   1;
            ];
        
        g_y_roots = [
            0.1   1;
            0.5   1;
            0.8   1;
            ];
        
        d_x_roots = [
            0.3   2;
            0.4   3;
            0.1   1;
            ];
        
        d_y_roots = [
            0.5   1;
            ];
        
        
       case 99
        f_x_roots = [
            2   1;
            3   1;
            4   1;
            9   1;
            ];
        f_y_roots = [
            2   1;
            3   1;
            4   1;
            5   1;
            6   1;
            ];
        
        g_x_roots = [
            2   1;
            4   1;
            7   1;
            9   1;
            ];
        g_y_roots = [
            2   1;
            5   1;
            7   1;
            ];
        
        d_x_roots = [
            2   1;
            4   1;
            9   1;
            ];
        
        d_y_roots = [
            2   1;
            5   1;
            ];
end

u_x_roots = getQuotient(f_x_roots,d_x_roots);
v_x_roots = getQuotient(g_x_roots,d_x_roots);
u_y_roots = getQuotient(f_y_roots,d_y_roots);
v_y_roots = getQuotient(g_y_roots,d_y_roots);




end

function u_x_roots = getQuotient(f_x_roots,d_x_roots)

num_roots_f_x = size(f_x_roots,1);
u_x_roots = [];
for i = 1:1:num_roots_f_x
    % get the root
    root = f_x_roots(i,1);
    % get multiplicity
    mult_f = f_x_roots(i,2);
    
    
    % Look up in d_x
    if ~isempty(find(d_x_roots(:,1) == root))
        
        [row_d,~] = find(d_x_roots(:,1) == root);
        % if root is found, get multiplicity
        mult_d = d_x_roots(row_d,2);
        
        % subtract to obtain multiplicity in u(x)
        mult_u = mult_f - mult_d;
        if mult_u > 0
            u_x_roots = [u_x_roots; root mult_u];
        end
    else
        u_x_roots = [u_x_roots; root mult_f];
    end
end
end