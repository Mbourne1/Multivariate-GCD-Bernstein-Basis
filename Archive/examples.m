function [f_x_roots,f_y_roots,...
    g_x_roots, g_y_roots,...
    u_x_roots, u_y_roots,...
    v_x_roots, v_y_roots,...
    d_x_roots, d_y_roots,...
    m1,m2,n1,n2,t1,t2] = examples(ex_num)




switch ex_num
    
    
    case 1
        % 1.1 get roots of polynomial f
        % % 1.1.1 get roots of polynomial f in terms of x
        f_x_roots = [
            5 1;
            1 1] ;
        % % 1.1.2 get roots of polynomial f in terms of y
        f_y_roots = [
            5 1];
        
        % 1.2 get roots of polynomial g
        % % 1.2.1 get roots of polynomial g in terms of x
        g_x_roots = [
            5 1;
            1 1];
        % % 1.2.2 get roots of polynomial g in terms of y
        g_y_roots = [-2 1];
        
        % 1.3 get roots of polynomial u
        % % 1.3.1 get roots of polynomial u in terms of x
        u_x_roots = [];
        % % 1.3.2 get roots of polynomial u in terms of y
        u_y_roots = [5 1];
        
        % 1.4 get roots of polynomial v
        % % 1.4.1 get roots of polyonmial v in terms of x
        v_x_roots = [
            ];
        % % 1.4.2 get roots of polynomial v in terms of y
        v_y_roots = [-2 1];
        
        % 1.5 get roots of polynomial d
        d_x_roots = [
            5 1;
            1 1];
        d_y_roots = [
            ];
        
    case 2
        
        % 1.1 get roots of polynomial f
        % % 1.1.1 get roots of polynomial f in terms of x
        f_x_roots = [
            0.1 1;
            0.2 1;
            0.3 1] ;
        % % 1.1.2 get roots of polynomial f in terms of y
        f_y_roots = [
            0.1 1;
            0.2 1;
            0.3 1];
        
        % 1.2 get roots of polynomial g
        % % 1.2.1 get roots of polynomial g in terms of x
        g_x_roots = [
            0.3 1;
            0.4 1;
            0.5 1];
        % % 1.2.2 get roots of polynomial g in terms of y
        g_y_roots = [
            0.3 1;
            0.4 1;
            0.5 1];
        
        % 1.3 get roots of polynomial u
        % % 1.3.1 get roots of polynomial u in terms of x
        u_x_roots = [
            0.1 1;
            0.2 1];
        
        % % 1.3.2 get roots of polynomial u in terms of y
        u_y_roots = [
            0.1 1;
            0.2 1];
        
        % 1.4 get roots of polynomial v
        % % 1.4.1 get roots of polyonmial v in terms of x
        v_x_roots = [
            0.4 1;
            0.5 1
            ];
        % % 1.4.2 get roots of polynomial v in terms of y
        v_y_roots = [
            0.4 1;
            0.5 1];
        
        % 1.5 get roots of polynomial d
        d_x_roots = [
            0.3 1
            ];
        d_y_roots = [
            0.3 1
            ];
        
        
    case 3
        % 1.1 get roots of polynomial f
        % % 1.1.1 get roots of polynomial f in terms of x
        f_x_roots = [
            5       1;
            2       1;
            1       1;
            ] ;
        % % 1.1.2 get roots of polynomial f in terms of y
        f_y_roots = [
            5       1
            1       1];
        
        % 1.2 get roots of polynomial g
        % % 1.2.1 get roots of polynomial g in terms of x
        g_x_roots = [
            5       1;
            2       1;
            1       1
            ];
        % % 1.2.2 get roots of polynomial g in terms of y
        g_y_roots = [
            2      1;
            1       1;
            ];
        
        % 1.3 get roots of polynomial u
        % % 1.3.1 get roots of polynomial u in terms of x
        u_x_roots = [];
        % % 1.3.2 get roots of polynomial u in terms of y
        u_y_roots = [
            5       1;
            ];
        
        % 1.4 get roots of polynomial v
        % % 1.4.1 get roots of polyonmial v in terms of x
        v_x_roots = [
            ];
        % % 1.4.2 get roots of polynomial v in terms of y
        v_y_roots = [
            2      1;
            ];
        
        % 1.5 get roots of polynomial d
        d_x_roots = [
            5       1;
            2       1;
            1       1;
            ];
        
        d_y_roots = [
            1   1;
            ];
        
    case 4
        % 1.1 get roots of polynomial f
        % % 1.1.1 get roots of polynomial f in terms of x
        f_x_roots = [
            0.5       1;
            0.3       1;
            0.15     1;
            0.1       1;
            ] ;
        % % 1.1.2 get roots of polynomial f in terms of y
        f_y_roots = [
            0.9       1;
            0.7       1;
            0.5       1;
            0.2       1;
            ];
        
        % 1.2 get roots of polynomial g
        % % 1.2.1 get roots of polynomial g in terms of x
        g_x_roots = [
            0.5       1;
            0.4       1;
            0.15      1;
            0.1       1;
            ];
        % % 1.2.2 get roots of polynomial g in terms of y
        g_y_roots = [
            0.9      1;
            0.7      1;
            0.3      1;
            0.2      1;
            ];
        
        % 1.3 get roots of polynomial u
        % % 1.3.1 get roots of polynomial u in terms of x
        u_x_roots = [
            0.3       1;
            ];
        % % 1.3.2 get roots of polynomial u in terms of y
        u_y_roots = [
            0.5       1;
            ];
        
        % 1.4 get roots of polynomial v
        % % 1.4.1 get roots of polyonmial v in terms of x
        v_x_roots = [
            0.4     1;
            ];
        % % 1.4.2 get roots of polynomial v in terms of y
        v_y_roots = [
            0.3     1;
            ];
        
        % 1.5 get roots of polynomial d
        d_x_roots = [
            0.5       1;
            0.15      1;
            0.1       1;
            ];
        d_y_roots = [
            0.9     1;
            0.7       1;
            0.2       1;
            ];
        
    case 5
        % 1.1 get roots of polynomial f
        % % 1.1.1 get roots of polynomial f in terms of x
        f_x_roots = [
            0.5       1;
            0.3       1;
            0.1       1;
            ] ;
        % % 1.1.2 get roots of polynomial f in terms of y
        f_y_roots = [
            0.3       1;
            0.7       1];
        
        % 1.2 get roots of polynomial g
        % % 1.2.1 get roots of polynomial g in terms of x
        g_x_roots = [
            0.5       1;
            0.1       1;
            ];
        % % 1.2.2 get roots of polynomial g in terms of y
        g_y_roots = [
            0.2      1;
            0.3       1;
            ];
        
        % 1.3 get roots of polynomial u
        % % 1.3.1 get roots of polynomial u in terms of x
        u_x_roots = [
            0.3       1;
            ];
        % % 1.3.2 get roots of polynomial u in terms of y
        u_y_roots = [
            0.7       1;
            ];
        
        % 1.4 get roots of polynomial v
        % % 1.4.1 get roots of polyonmial v in terms of x
        v_x_roots = [
            ];
        % % 1.4.2 get roots of polynomial v in terms of y
        v_y_roots = [
            0.2 1
            ];
        
        % 1.5 get roots of polynomial d
        d_x_roots = [
            0.5       1;
            0.1       1];
        
        d_y_roots = [
            0.3       1
            ];
        
        
    case 6
        % 1.1 get roots of polynomial f
        % % 1.1.1 get roots of polynomial f in terms of x
        f_x_roots = [
            0.5       1;
            0.3       1;
            0.1       1;
            0.7       1;
            ] ;
        % % 1.1.2 get roots of polynomial f in terms of y
        f_y_roots = [
            0.3       1;
            0.5       1];
        
        % 1.2 get roots of polynomial g
        % % 1.2.1 get roots of polynomial g in terms of x
        g_x_roots = [
            0.5       1;
            0.1       1;
            ];
        % % 1.2.2 get roots of polynomial g in terms of y
        g_y_roots = [
            -0.2      1;
            0.3       1;
            ];
        
        % 1.3 get roots of polynomial u
        % % 1.3.1 get roots of polynomial u in terms of x
        u_x_roots = [
            0.3       1;
            0.7       1;
            ];
        % % 1.3.2 get roots of polynomial u in terms of y
        u_y_roots = [
            0.5       1;
            ];
        
        % 1.4 get roots of polynomial v
        % % 1.4.1 get roots of polyonmial v in terms of x
        v_x_roots = [
            ];
        % % 1.4.2 get roots of polynomial v in terms of y
        v_y_roots = [
            -0.2 1
            ];
        
        % 1.5 get roots of polynomial d
        d_x_roots = [
            0.5       1;
            0.1       1];
        
        d_y_roots = [
            0.3       1
            ];
        
     case 8
        % 1.1 get roots of polynomial f
        % % 1.1.1 get roots of polynomial f in terms of x
        f_x_roots = [
            0.5       1;
            0.3       1;
            0.15     1;
            0.1       1;
            0.8     1;
            ] ;
        % % 1.1.2 get roots of polynomial f in terms of y
        f_y_roots = [
            0.9       1;
            0.8       1;
            0.7       1;
            0.5       1;
            0.2       1;
            ];
        
        % 1.2 get roots of polynomial g
        % % 1.2.1 get roots of polynomial g in terms of x
        g_x_roots = [
            0.6     1;
            0.5       1;
            0.4       1;
            0.15      1;
            0.1       1;
            ];
        % % 1.2.2 get roots of polynomial g in terms of y
        g_y_roots = [
            0.9      1;
            0.7      1;
            0.3      1;
            0.2      1;
            0.1     1;
            ];
        
        % 1.3 get roots of polynomial u
        % % 1.3.1 get roots of polynomial u in terms of x
        u_x_roots = [
            0.3       1;
            0.8     1;
            ];
        % % 1.3.2 get roots of polynomial u in terms of y
        u_y_roots = [
            0.5       1;
            0.8     1;
            ];
        
        % 1.4 get roots of polynomial v
        % % 1.4.1 get roots of polyonmial v in terms of x
        v_x_roots = [
            0.4     1;
            0.6     1;
            ];
        % % 1.4.2 get roots of polynomial v in terms of y
        v_y_roots = [
            0.3     1;
            0.1     1;
            ];
        
        % 1.5 get roots of polynomial d
        d_x_roots = [
            0.5       1;
            0.15      1;
            0.1       1;
            ];
        d_y_roots = [
            0.9     1;
            0.7       1;
            0.2       1;
            ];
end



% Get the degree of m1 = number of roots of f in terms of x
if isempty(f_x_roots)
    m1 = 0;
else
    m1 = sum(f_x_roots(:,2));
end
% Get the degree of m2 = number of roots of f in terms of y
if isempty(f_y_roots)
    m2 = 0;
else
    m2 = sum(f_y_roots(:,2));
end

% Get the degree of n1 = number of roots of g in terms of x
if isempty(g_x_roots)
    n1 = 0;
else
    n1 = sum(g_x_roots(:,2));
end

% Get the degree of n2 = number of roots of g in terms of y
if isempty(g_y_roots)
    n2 = 0;
else
    n2 = sum(g_y_roots(:,2));
end

% Get the degree of t1 = number of roots of d in terms of y
if isempty(d_x_roots)
    t1 = 0;
else
    t1 = sum(d_x_roots(:,2));
end

% Get the degree of t2 = number of roots of d in terms of y
if isempty(d_y_roots)
    t2 = 0;
else
    t2 = sum(d_y_roots(:,2));
end
end
