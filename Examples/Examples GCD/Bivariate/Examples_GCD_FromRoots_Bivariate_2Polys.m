function [fxy, gxy, dxy, uxy, vxy, m, n, t] = Examples_GCD_FromRoots_Bivariate_2Polys(ex_num)
% Get a GCD example given an example number
%
% % Inputs.
%
% ex_num : (String) Example Number as a string
%
% % Outputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y) 
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% dxy : (Matrix) Coefficients of polynomial d(x,y), the greatest common
% divisor of f(x,y) and g(x,y)
%
% uxy : (Matrix) Coefficients of polynomial u(x,y), the quotient polynomial
% such that f(x,y)/u(x,y) = d(x,y) 
%
% vxy : (Matrix) Coefficients of polynomial v(x,y), the quotient polynomial
% such that g(x,y)/v(x,y) = d(x,y)
%
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% t : (Int) Total degree of d(x,y)

f_roots_x = [];
f_roots_y = [];
f_roots_xy = [];

g_roots_x = [];
g_roots_y = [];
g_roots_xy = [];

d_roots_x = [];
d_roots_y = [];
d_roots_xy = [];

u_roots_x = [];
u_roots_y = [];
u_roots_xy = [];

v_roots_x = [];
v_roots_y = [];
v_roots_xy = [];



switch ex_num
    
    case 'Custom'
        
        prompt = 'Enter the degree of Polynomial f(x) :';
        m = input(prompt);
        
        prompt = 'Enter the degree of Polynomial g(x) :';
        n = input(prompt);
        
        prompt = 'Enter the degree of Polynomial d(x) :';
        t = input(prompt);
        
        intvl_low = -1;
        intvl_high = 1;
        
        [f_roots_x,g_roots_x] = BuildRandomPolynomials(m,n,t,intvl_low, intvl_high);
        [f_roots_y,g_roots_y] = BuildRandomPolynomials(m,n,t,intvl_low, intvl_high);
        
        d_roots_x = GetDivisor(f_roots_x,g_roots_x);
        d_roots_y = GetDivisor(f_roots_y,g_roots_y);
        
        u_roots_x = GetQuotient(f_roots_x,d_roots_x);
        u_roots_y = GetQuotient(f_roots_y,d_roots_y);
        
        v_roots_x = GetQuotient(g_roots_x,d_roots_x);
        v_roots_y = GetQuotient(g_roots_y,d_roots_y);
        
    case 'template'
        
        % Roots of Polynomial f(x,y)
        f_roots_x = [...
            ];
        f_roots_y = [...
            ];
        f_roots_xy = [];
        
        % Roots of polynomial g(x,y)
        g_roots_x = [...
            ];
        g_roots_y = [...
            ];
        g_roots_xy = [];
        
        % Roots of polynomial d(x,y)
        d_roots_x = [...
            ];
        d_roots_y = [...
            ];
        d_roots_xy = [];
        
        % Roots of polynomial u(x,y)
        u_roots_x = [...
            ];
        u_roots_y = [...
            ];
        u_roots_xy = [];
        
        
        % Roots of polynomial v(x,y)
        v_roots_x = [...
            ];
        v_roots_y = [...
            ];
        v_roots_xy = [];
        
        m = 100;
        n = 100;
        t = 100;
        
        
        
    case '1'
        % Example 1
        % From the file Bivariate - Bernstein Basis - GCD Examples.tex
        % This is a separable example
        
        f_roots_x = [...
            0.7 1
            0.5 1
            0.2 1
            ];
        f_roots_y = [...
            -0.2    1
            ];
        
        % Root of g with respect to x
        g_roots_x = ...
            [
            0.7 1
            ];
        
        % Roots of g with respect to y
        g_roots_y = ...
            [
            -0.2    1
            -0.9    1
            ];
        
        % Roots of u with respect to x
        u_roots_x = ...
            [
            0.5 1
            0.2 1
            ];
        
        % Roots of u with respect to y
        u_roots_y = ...
            [
            ];
        
        % Roots of v with respect to x
        v_roots_x = ...
            [
            ];
        % Roots of v with respect to y
        v_roots_y = ...
            [
            -0.9    1
            ];
        
        % Roots of d with respcet to x
        d_roots_x = ...
            [
            0.7 1;
            ];
        % Roots of y with respect to y
        d_roots_y = ...
            [
            -0.2    1
            ];
        
        
        m = 4;
        n = 3;
        t = 2;
        
    case '2'
        %
        % Example 2
        % From the file Bivariate - Bernstein Basis - GCD Examples.tex
        % This is a non separable example
        %
        
        % Roots of f with respect to x
        f_roots_x = [...
            0.7     1
            0.5     1
            0.2     1
            ];
        
        % Roots of f with respect to y
        f_roots_y = [...
            -0.2    1;
            ];
        
        % Roots of f with respect to xy
        f_roots_xy{1,1} = ...
            [
            -0.5    0.5
            0.5     1.5
            ];
        
        % Roots of g with respect to x
        g_roots_x = [...
            0.7 1
            ];
        
        % Roots of g with respect to y
        g_roots_y = [...
            -0.2    1
            -0.9    1
            ];
        
        % Roots of g with respect to xy)
        g_roots_xy{1,1} = ...
            [
            -0.5    0.5
            0.5     1.5
            ];
        
        % Roots of u with respect to x
        u_roots_x = ...
            [
            0.5 1
            0.2 1
            ];
        
        % Roots of v with respect to x
        v_roots_x = ...
            [
            ];
        
        % Roots of v with respect to y
        v_roots_y  = ...
            [
            -0.9    1
            ];
        
        % Roots of d with respect to x
        d_roots_x = ...
            [
            0.7     1
            ];
        % Roots of d with respect to y
        d_roots_y = ...
            [
            -0.2    1
            ];
        
        % Roots of d with respect to xy
        d_roots_xy{1,1} = ...
            [
            -0.5    0.5
            0.5     1.5
            ];
        
        % Degrees
        m = 5;
        n = 4;
        t = 3;
        
    case '3'
        %
        % Example 3
        % From the file Bivariate - Bernstein Basis - GCD Examples.tex
        % This is a non separable example
        % An extra root in f_x
        
        % Roots of f with respect to x
        f_roots_x =...
            [
            0.7     1;
            0.5     1;
            0.2     1;
            0.1    1;
            ];
        % Roots of f with respect to y
        f_roots_y = ...
            [
            0.2     1;
            0.3     1;
            ];
        
        
        % non separable root (x+y-0.5)
        f_roots_xy{1,1} = ...
            [
            -0.5    0.5
            0.5     1.5
            ];
        
        % Roots of polynomial g(x,y)
        % Roots of g with respect to x
        g_roots_x = ...
            [
            0.7     1;
            0.9     1;
            ];
        % Roots of g with respect to y
        g_roots_y = ...
            [
            0.2    1
            0.9    1
            ];
        
        % Roots of g with respect to x and y
        g_roots_xy{1,1} = ...
            [
            -0.5    0.5
            0.5     1.5
            ];
        
        % Roots of polynomial u(x,y)
        % Roots of u with respect to x
        u_roots_x =...
            [
            0.5     1
            0.2     1
            0.1     1
            ];
        % Roots of u with respect to y
        u_roots_y = ...
            [
            0.3     1
            ];
        
        % Roots of Polynomial v(x,y)
        % Roots of v with respect to x
        v_roots_x = ...
            [
            0.9     1;
            ];
        % Roots of v with respect to y
        v_roots_y = ...
            [
            0.9    1;
            ];
        
        % Roots of d with respect to x
        d_roots_x = ...
            [
            0.7     1
            ];
        % Roots of d with respect to y
        d_roots_y = ...
            [
            0.2    1
            ];
        % Roots of d with respect to xy
        d_roots_xy{1,1} = ...
            [
            -0.5    0.5
            0.5     1.5
            ];
        
        
        m = 7;
        n = 5;
        t = 3;
        
        
    case '4'
        
        % Roots of f with respect to x
        f_roots_x = ...
            [
            2.4   1
            1.5   2
            1.9   3
            ];
        % Roots of f with respect to y
        f_roots_y = ...
            [
            0.3   1
            1.1   2
            1.7   2
            ];
        % Roots of f with respect to xy
        f_roots_xy{1,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        
        f_roots_xy{2,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        f_roots_xy{3,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        % Roots of g with respect to x
        g_roots_x = ...
            [
            2.4   1
            1.5   2
            2.9   3
            ];
        
        % Roots of g with respect to y
        g_roots_y = ...
            [
            1.9   1
            1.1   2
            ];
        
        % Roots of f with respect to xy
        g_roots_xy{1,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        g_roots_xy{2,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        g_roots_xy{3,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        g_roots_xy{4,1} =...
            [
            1.7     2.7;
            2.7     3.7;
            ];
        
        
        % Roots of d with respect to x
        d_roots_x = ...
            [
            2.4   1
            1.5   2
            ];
        % Roots of d with respect to y
        d_roots_y = ...
            [
            1.1   2;
            ];
        
        % Roots of d with respect to xy
        d_roots_xy{1,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        d_roots_xy{2,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        d_roots_xy{3,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        
        
        % Roots of u with respect to x
        u_roots_x = ...
            [
            1.9   3
            ];
        
        % Roots of u with respect to y
        u_roots_y = ...
            [
            0.3   1;
            1.7   2;
            ];
        
        
        % Roots of v with respect to x
        v_roots_x = ...
            [
            2.9   3
            ];
        
        % Roots of v with respect to y
        v_roots_y = ...
            [
            1.9   1
            ];
        
        % Roots of v with respect to xy
        v_roots_xy = ...
            [
            1.7     2.7;
            2.7     3.7;
            ];
        
        
        m = 17;
        n = 17;
        t = 11;
        
        
    case '5'
        %% EXAMPLE :
        % Roots of f with respect to x
        f_roots_x = ...
            [
            1.4   1
            0.5   2
            0.9   3
            ];
        % Roots of f with respect to y
        f_roots_y = ...
            [
            1.3   1
            1.1   2
            1.7   2
            ];
        % Roots of f with respect to xy
        f_roots_xy{1,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        f_roots_xy{2,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        f_roots_xy{3,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        
        % Roots of g with respect to x
        g_roots_x = ...
            [
            1.4   1
            0.5   2
            0.4   3
            ];
        
        % Roots of g with respect to y
        g_roots_y = ...
            [
            1.9   1
            1.1   2
            ];
        
        % Roots of f with respect to xy
        g_roots_xy{1,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        g_roots_xy{2,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        g_roots_xy{3,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        g_roots_xy{4,1} =...
            [
            1.7     2.7;
            2.7     3.7;
            ];
        
        
        % Roots of d with respect to x
        d_roots_x = ...
            [
            1.4   1
            0.5   2
            ];
        % Roots of d with respect to y
        d_roots_y = ...
            [
            1.1   2;
            ];
        
        % Roots of d with respect to xy
        d_roots_xy{1,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        d_roots_xy{2,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        d_roots_xy{2,1} = ...
            [
            0.1     1.1;
            1.1     2.1
            ];
        
        
        % Roots of u with respect to x
        u_roots_x = ...
            [
            0.9   3
            ];
        
        % Roots of u with respect to y
        u_roots_y = ...
            [
            1.3   1;
            1.7   2;
            ];
        
        
        % Roots of v with respect to x
        v_roots_x = ...
            [
            0.4   3
            ];
        
        % Roots of v with respect to y
        v_roots_y = ...
            [
            1.9   1
            ];
        
        % Roots of v with respect to xy
        v_roots_xy = ...
            [
            1.7     2.7;
            2.7     3.7;
            ];
        
        
        m = 14;
        n = 13;
        t = 8;
    case '6'
        % Example with factorisation
        % (x+y+0.1) (x-0.2) (y-0.3)
        
        % Roots of f with respect to x
        f_roots_x = ...
            [
            0.2   1;
            ];
        
        % Roots of f with respect to y
        f_roots_y = ...
            [
            0.3   1;
            ];
        
        % Roots of f with respect to x and y
        f_roots_xy{1,1} = [...
            0.1  1.1;
            1.1  2.1
            ];
        
        
        % Roots of g with respect to x
        g_roots_x = ...
            [
            ];
        % Roots of g with respect to y
        g_roots_y = ...
            [
            0.3     1;
            ];
        
        % Roots in both x and y
        g_roots_xy{1,1} = ...
            [
            0.3     1.3
            1.3     2.3
            ];
        
        % Roots of d with respect to x
        d_roots_x = [];
        
        % Roots of d with respect to y
        d_roots_y = ...
            [
            0.3     1;
            ];
        
        % Roots of u with respect to x
        u_roots_x = ...
            [
            0.2   1;
            ];
        
        % Roots of u with respect to xy
        u_roots_xy{1,1} = [...
            0.1  1.1;
            1.1  2.1
            ];
        
        
        
        % Roots of v(x,y)
        % Roots in both x and y
        v_roots_xy{1,1} = [...
            0.3     1.3
            1.3     2.3
            ];
        
        
        m = 4;
        n = 3;
        t = 1;
        %
    case '7'
        % Example 18 in Power File
        
        % Roots of f with respect to x
        f_roots_x = ...
            [
            0.573429421     1;
            2.2175          3;
            0.8765          2;
            0.1759262       1;
            1.1057853       1;
            ];
        
        % Roots of f with respect to xy
        f_roots_xy{1,1} = [...
            0.1  1.1;
            1.1  2.1
            ];
        
        % Roots of g with respect to x
        g_roots_x = ...
            [
            0.573429421     1;
            0.7891          1;
            0.1234          1;
            0.8765          2;
            0.0175267       1;
            ];
        
        % Roots of g with respect to xy
        g_roots_xy{1,1} = [...
            0.1  1.1;
            1.1  2.1
            ];
        
        % Roots of d with respect to x
        d_roots_x = ...
            [
            0.573429421     1;
            0.8765          2;
            ];
        
        % Roots of f with respect to xy
        d_roots_xy{1,1} = [...
            0.1  1.1;
            1.1  2.1
            ];
        
        % Roots of u with respect to x
        u_roots_x = ...
            [
            2.2175          3;
            0.1759262       1;
            1.1057853       1;
            ];
        
        % Roots of v with respect to x
        v_roots_x = ...
            [
            0.7891          1;
            0.1234          1;
            0.0175267       1;
            ];
        
        
        m = 9;
        n = 7;
        t = 4;
        
    case '8'
        % Examle 19 in power file
        
        % Roots of f with respect to x
        f_roots_x = ...
            [
            0.573429421 1;
            2.2175      3;
            0.8765      2;
            0.1759262   1;
            1.1057853   1;
            ];
        % Roots of f with respect to y
        f_roots_y = ...
            [
            0.424242    3
            ];
        % Roots of g with respect to x
        g_roots_x = ...
            [
            0.573429421 1;
            0.7891      1;
            0.1234      1;
            0.8765      2;
            0.0175267   1;
            ];
        
        % Roots of g with respect to y
        g_roots_y = ...
            [
            0.424242    3;
            ];
        
        %Roots of d with respect to x
        d_roots_x = ...
            [
            0.573429421 1;
            0.8765      2;
            ];
        % Roots of d with respect to y
        d_roots_y = ...
            [
            0.424242    3
            ];
        
        % Roots of u with respect to x
        u_roots_x = ...
            [
            2.2175      3;
            0.1759262   1;
            1.1057853   1;
            ];
        
        % Roots of v with respect to x
        v_roots_x = ...
            [
            0.7891      1;
            0.1234      1;
            0.0175267   1;
            ];
        
        m = 11;
        n = 9;
        t = 8;
        
    case '9'
        % Roots of f with respect to x
        f_roots_x = ...
            [
            0.573429421 1;
            2.2175      3;
            0.8765      4;
            0.1759262   1;
            1.1057853   1;
            ];
        
        % Roots in f with respect to y
        f_roots_y = ...
            [
            0.424242    3;
            ];
        % Roots in g with respect to x
        g_roots_x = ...
            [
            0.573429421 1
            0.7891      1
            0.1234      1
            0.8765      4
            0.0175267   1
            ];
        % Roots in g with respect to y
        g_roots_y = ...
            [
            0.424242    3
            ];
        
        % Roots in d with respect to x
        d_roots_x = ...
            [
            0.573429421 1;
            0.8765      4;
            ];
        % Roots in d with respect to y
        d_roots_y = ...
            [
            0.424242    3;
            ];
        
        % Roots in u with respect to x
        u_roots_x = ...
            [
            2.2175      3;
            0.1759262   1;
            1.1057853   1;
            ];
        % Roots in v with respect to x
        v_roots_x = ...
            [
            0.7891      1;
            0.1234      1;
            0.0175267   1;
            ];
        
        m = 13;
        n = 11;
        t = 8;
        % Corresponds to case 20 in power basis examples
    case '10'
        
        % Roots of f with respect to x
        f_roots_x = ...
            [
            0.573429421 1
            2.2175      3
            0.8765      4
            0.1759262   1
            1.1057853   4
            ];
        
        % Roots of f with respect to y
        f_roots_y = ...
            [
            0.424242    3;
            ];
        % Roots of g with respect to x
        g_roots_x = ...
            [
            0.573429421 1
            0.7891      1
            0.1234      1
            0.8765      4
            0.0175267   1
            ];
        
        % Roots of g with respect to y
        g_roots_y = ...
            [
            0.424242    3;
            ];
        % Roots of d with respect to x
        d_roots_x = ...
            [
            0.573429421 1;
            0.8765      4;
            ];
        
        % Roots of g with respect to y
        d_roots_y = ...
            [
            0.424242    3;
            ];
        % Roots of u with respect to x
        u_roots_x = ...
            [
            2.2175  3
            0.1759262   1
            1.1057853   4
            ];
        % roots of v with respect to x
        v_roots_x = ...
            [
            0.7891      1;
            0.1234      1;
            0.0175267   1;
            ];
        
        m = 16;
        n = 11;
        t = 8;
        
        
    case '11'
        % From winkler - Two methods...
        
        f_roots_x = ...
            [
            0.5161 5;
            -7.1052 5;
            0.1132 3
            ];
        
        g_roots_x = ...
            [
            0.5161 5;
            -7.1052 5;
            2.0476  7;
            -8.8614 7;
            ];
        
        d_roots_x = ...
            [
            0.5161  5
            -7.1052  5
            ];
        
        u_roots_x = ...
            [
            0.1132  3
            ];
        
        v_roots_x = ...
            [
            2.0476  7;
            -8.8614  7;
            ];
        
        m = 13;
        n = 24;
        t = 10;
        
        
    case '12'
        % From winkler - Two methods...
        % With added roots in y
        
        % Roots of f with respect to x
        f_roots_x =...
            [
            -0.5161 5;
            0.1052 5;
            -0.1132 3;
            ];
        
        % Roots of f with respect to y
        f_roots_y = ...
            [
            0.7725  3;
            ];
        
        % Roots of g with respect to x
        g_roots_x = ...
            [
            -0.5161     5;
            0.1052     5;
            2.0476     7;
            0.8614     7;
            ];
        
        % Roots of g with respect to y
        g_roots_y = ...
            [
            0.7725  3;
            ];
        
        % Roots of d with respect to x
        d_roots_x = ...
            [
            -0.5161  5;
            0.1052  5;
            ];
        
        % Roots of d in y
        d_roots_y = ...
            [
            0.7725  3;
            ];
        
        % Roots of u in x
        u_roots_x = ...
            [
            -0.1132  3;
            ];
        
        % Roots of v in x
        v_roots_x = ...
            [
            2.0476  7
            0.8614  7
            ];
        
        m = 16;
        n = 27;
        t = 13;
    otherwise
        
        
        error('not a valid example number')
        
        
end

f_roots_x = mult_roots_x(f_roots_x);
f_roots_y = mult_roots_y(f_roots_y);

g_roots_x = mult_roots_x(g_roots_x);
g_roots_y = mult_roots_y(g_roots_y);

d_roots_x = mult_roots_x(d_roots_x);
d_roots_y = mult_roots_y(d_roots_y);

u_roots_x = mult_roots_x(u_roots_x);
u_roots_y = mult_roots_y(u_roots_y);

v_roots_x = mult_roots_x(v_roots_x);
v_roots_y = mult_roots_y(v_roots_y);

f_roots = [f_roots_x ; f_roots_y ; f_roots_xy];
g_roots = [g_roots_x ; g_roots_y ; g_roots_xy];
d_roots = [d_roots_x ; d_roots_y ; d_roots_xy];
u_roots = [u_roots_x ; u_roots_y ; u_roots_xy];
v_roots = [v_roots_x ; v_roots_y ; v_roots_xy];




fxy = BuildPoly_NonSeparable(f_roots);
gxy = BuildPoly_NonSeparable(g_roots);
uxy = BuildPoly_NonSeparable(u_roots);
vxy = BuildPoly_NonSeparable(v_roots);
dxy = BuildPoly_NonSeparable(d_roots);

[m1, m2] = GetDegree_Bivariate(fxy);

[n1, n2] = GetDegree_Bivariate(gxy);

[t1, t2] = GetDegree_Bivariate(dxy);


end



function cellArr = mult_roots_x(root_mult_mat)
% Given the root and multiplicity matrix return 

%                    _______________
% root_mult_mat =   | root  |  mult |
%                   |       |       |
%                   |_______|_______|


count = 1;

cellArr = {};

% for each root in the array
[num_roots,~ ] = size(root_mult_mat);

for i = 1:1:num_roots
    % get the multiplicity of the root
    mult = root_mult_mat(i,2);
    % get the root
    root = root_mult_mat(i,1);
    for j = 1:1:mult
        cellArr{count,1} = root_x(root);
        count = count + 1;
    end
end

end

function cellArr = mult_roots_y(root_mult_mat)
% given the root and multiplicity matrix

%                    _______________
% root_mult_mat =   | root  |  mult |
%                   |       |       |
%                   |_______|_______|


count = 1;
cellArr = {};
% for each root in the array
[num_roots,~ ] = size(root_mult_mat);

for i = 1:1:num_roots
    % get the multiplicity of the root
    mult = root_mult_mat(i,2);
    % get the root;
    root = root_mult_mat(i,1);
    for j = 1:1:mult
        cellArr{count,1} = root_y(root);
        count = count + 1;
    end
end
end