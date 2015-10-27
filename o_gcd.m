function [dxy_matrix] = o_gcd()
% return the GCD of two polynomials

%%
global bool_preproc
bool_preproc = input('Do you wish to apply preprocessing operations (y)es/(n)o :  ','s');
global bool_noise
bool_noise = input('Do you wish to add noise to the coefficients (y)es/(n)o :  ','s');

global bool_Q
bool_Q = input('Do you wish to include Q in the coefficient matrix (y)es/(n)o :  ','s');


%% Get Example details

ex_type = input('Example Type (s)eparable or (ns) non separable : ','s');

switch ex_type
    case 's'
                
        [fxy_matrix,gxy_matrix,...
            uxy_matrix,vxy_matrix,...
            dxy_matrix,m,n,t,t1,t2] = Examples(ex_num);
               
        fprintf('Exact Degree of GCD : %i \n', d)
        fprintf('Exact t1 : %i \n',t1)
        fprintf('Exact t2 : %i \n',t2)
        
        
    case 'ns'
        
        
        example = input('Enter Example Number: ', 's');
        [fxy_mtrx_exct, gxy_mtrx_exct,...
            uxy_mtrx_exct,vxy_mtrx_exct,...
            dxy_mtrx_exct,...
            m,m1,m2,...
            n,n1,n2,...
            t_exct,t1_exct,t2_exct] = Examples_NonSeparableWithSolutions(example);
        
        
        fprintf('\n')
        fprintf('----------------------------------------------------------------\n')
        fprintf('Input Polynomials Degrees:\n')
        fprintf('m  : %i \n',m)
        fprintf('m1 : %i \n',m1)
        fprintf('m2 : %i \n\n',m2)
        fprintf('n  : %i \n',n)
        fprintf('n1 : %i \n',n1)
        fprintf('n2 : %i \n\n',n2)
        
        fprintf('t  : %i \n',t_exct)
        fprintf('t1 : %i \n',t1_exct)
        fprintf('t2 : %i \n',t2_exct)
        fprintf('----------------------------------------------------------------\n')
        fprintf('\n')
        fprintf('----------------------------------------------------------------\n')
        fprintf('\n')
        fprintf('Polynomial d(x,y)\n')
        dxy_mtrx_exct
        fprintf('----------------------------------------------------------------\n')
        
        
        
end

%%
% Take a working copy of f and g
fxy_matrix_working = fxy_mtrx_exct;
gxy_matrix_working = gxy_mtrx_exct;

%% Noise
switch bool_noise
    case 'y'
        % Add noise to the coefficients of f and g
        
        % Set minimum and maximum noise levels
        emin_input = input('Enter signal:noise ratio : ','s');
        emax_input = input('Enter singal:noise ratio : ','s');
        
        emin = str2num(emin_input);
        emax = str2num(emax_input);
        
        % Add noise to the coefficients of f and g
        fxy_matrix_working = Noise(fxy_matrix_working,emin,emax);
        gxy_matrix_working = Noise(gxy_matrix_working,emin,emax);
        
        fprintf('\n')
        fprintf('Noise added to coefficients of fxy\n')
        noise_matrix = fxy_matrix_working - fxy_mtrx_exct
        fprintf('\n')
    case 'n'
        
end


%% given the two polynomials fxy and gxy - plot them
ans = input('Do you wish to plot fxy and gxy? (y/n) :  ','s');
if ans == 'y'
    plot_fxy(fxy_matrix_working);
    plot_fxy(gxy_matrix_working);
else
    % Dont plot
end

[uxy_matrix, vxy_matrix, dxy_matrix] = o1(fxy_matrix_working,gxy_matrix_working,...
    m,n,...
    t_exct,t1_exct,t2_exct,...
    uxy_mtrx_exct,vxy_mtrx_exct,dxy_mtrx_exct,...
    fxy_mtrx_exct,gxy_mtrx_exct)


end