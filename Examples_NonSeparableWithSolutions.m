function [fxy_mtrx_exct,gxy_mtrx_exct,...
    uxy_mtrx_exct, vxy_mtrx_exct,...
    dxy_mtrx_exct,...
    m,m1,m2,...
    n,n1,n2,...
    t,t1,t2] = Examples_NonSeparableWithSolutions(ex)

% NOTE -THESE EXAMPLES ASSUME SMALLEST POWER FIRST


switch ex
    case '0'
        % Separable example
        % Bernstein conversion of example 0 from pwr basis examples
        fxy_mtrx_exct = [...
            12 12 9 6;
            8/3 8/3 2 6;
            0 0 0 6;
            0 0 0 6
            ]
        gxy_mtrx_exct = [...
            -28 -28 -21 0;
            -8 -8 -6 0;
            0 0 0 0
            ];
        dxy_mtrx_exct = [...
            4 4 3 0;
            0 0 0 0
            ];
        uxy_mtrx_exct = [...
            3;
            1;
            0
            ];
        vxy_mtrx_exct = [...
            -7
            -4
            ];
        
        m = 6;
        n = 5;
        t = 4;
    case '1'
        fxy_mtrx_exct = [...
            1   2   3   0;
            0   0   0   0
            ]
        gxy_mtrx_exct = [...
            1   2   3   0;
            ];
        dxy_mtrx_exct = [...
            1   2   3   0
            ];
        uxy_mtrx_exct = [...
            1;
            0
            ];
        vxy_mtrx_exct = [...
            1
            ];
        
        m = 4;
        n = 3;
        t = 3;
        
        
        
        
    case '20'
        % Example 1
        % From the file Bivariate - Bernstein Basis - GCD Examples.tex
        % This is a separable example
        
        
        roots_f{1,1} = ...
            [
            -0.7;
            0.3
            ];
        roots_f{2,1} = ...
            [
            -0.5;
            0.5
            ];
        roots_f{3,1} = ...
            [
            -0.2
            0.8
            ];
        % % Roots in f with respect to y
        roots_f{4,1} = ...
            [
            0.2    1.2
            ];
        
        roots_g{1,1} = ...
            [
            -0.7;
            0.3
            ];
        roots_g{2,1} = ...
            [
            0.2  1.2
            ];
        roots_g{3,1} = ...
            [
            0.9  1.9
            ];
        
        roots_u{1,1} = ...
            [
            -0.5;
            0.5
            ];
        roots_u{2,1} = ...
            [
            -0.2
            0.8
            ];
        roots_v{1,1} = ...
            [
            0.9  1.9
            ];
        roots_d{1,1} = ...
            [
            -0.7;
            0.3
            ];
        roots_d{2,1} = ...
            [
            0.2 1.2
            ];
        
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f)
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g)
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u)
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v)
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d)
        
        m = 3;
        n = 3;
        t = 2;
        
    case '21'
        %%
        % Example 2
        % From the file Bivariate - Bernstein Basis - GCD Examples.tex
        % This is a non separable example
        %%
        %%
        roots_f{1,1} = ...
            [
            -0.7;
            0.3
            ];
        roots_f{2,1} = ...
            [
            -0.5;
            0.5
            ];
        roots_f{3,1} = ...
            [
            -0.2
            0.8
            ];
        % non separable root (x+y-0.5)
        roots_f{4,1} = ...
            [
            -0.5    0.5
            0.5     1.5
            ];
        % % Roots in f with respect to y
        roots_f{5,1} = ...
            [
            0.2    1.2
            ];
        
        roots_g{1,1} = ...
            [
            -0.7;
            0.3
            ];
        roots_g{2,1} = ...
            [
            0.2  1.2
            ];
        roots_g{3,1} = ...
            [
            0.9  1.9
            ];
        % Non separable root (x+y-0.5)
        roots_g{4,1} = ...
            [
            -0.5    0.5
            0.5     1.5
            ];
        
        roots_u{1,1} = ...
            [
            -0.5;
            0.5
            ];
        roots_u{2,1} = ...
            [
            -0.2
            0.8
            ];
        roots_v{1,1} = ...
            [
            0.9  1.9
            ];
        roots_d{1,1} = ...
            [
            -0.7;
            0.3
            ];
        roots_d{2,1} = ...
            [
            0.2 1.2
            ];
        roots_d{3,1} = ...
            [
            -0.5    0.5
            0.5     1.5
            ];
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f)
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g)
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u)
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v)
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d)
        
        m = 4;
        n = 4;
        t = 3;
        
        
        case '22'
        %%
        % Example 3
        % From the file Bivariate - Bernstein Basis - GCD Examples.tex
        % This is a non separable example
        % An extra root in f_x
        %%
        %%
        roots_f{1,1} = ...
            [
            -0.7;
            0.3
            ];
        roots_f{2,1} = ...
            [
            -0.5;
            0.5
            ];
        roots_f{3,1} = ...
            [
            -0.2
            0.8
            ];
        % non separable root (x+y-0.5)
        roots_f{4,1} = ...
            [
            -0.5    0.5
            0.5     1.5
            ];
        % % Roots in f with respect to y
        roots_f{5,1} = ...
            [
            0.2    1.2
            ];
        
        roots_f{6,1} = ...
            [
            0.1;
            1.1
            ];
        roots_f{7,1} = ...
            [
            0.75654;
            1.75654
            ];
        roots_f{8,1} = ...
            [
            0.3759  1.3759
            ];
        roots_f{9,1} = ...
            [
            0.56789; 
            1.56789
            ]
        %%
        roots_g{1,1} = ...
            [
            -0.7;
            0.3
            ];
        roots_g{2,1} = ...
            [
            0.2  1.2
            ];
        roots_g{3,1} = ...
            [
            0.9  1.9
            ];
        % Non separable root (x+y-0.5)
        roots_g{4,1} = ...
            [
            -0.5    0.5
            0.5     1.5
            ];
        roots_g{5,1} = ...
            [
            -0.9;
            0.1
            ];
        
        %%
        roots_u{1,1} = ...
            [
            -0.5;
            0.5
            ];
        roots_u{2,1} = ...
            [
            -0.2
            0.8
            ];
        roots_u{3,1} = ...
            [
            0.1;
            1.1
            ];
        roots_u{4,1} = ...
            [
            0.75654;
            1.75654
            ];
        roots_u{5,1} = ...
            [
            0.3759  1.3759
            ];
        roots_u{6,1} = ...
            [
            0.56789; 
            1.56789
            ]
        %%
        roots_v{1,1} = ...
            [
            0.9  1.9
            ];
        roots_v{2,1} = ...
            [
            -0.9;
            0.1
            ];
        
        
        %%
        roots_d{1,1} = ...
            [
            -0.7;
            0.3
            ];
        roots_d{2,1} = ...
            [
            0.2 1.2
            ];
        roots_d{3,1} = ...
            [
            -0.5    0.5
            0.5     1.5
            ];
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f)
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g)
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u)
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v)
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d)
        
        m = 7;
        n = 5;
        t = 3;
        
        
     
        
end

[r,c] = size(fxy_mtrx_exct);
m1 = r -1;
m2 = c -1;

[r,c] = size(gxy_mtrx_exct);
n1 = r -1;
n2 = c -1;

[r,c] = size(dxy_mtrx_exct);
t1 = r-1;
t2 = c-1;


end