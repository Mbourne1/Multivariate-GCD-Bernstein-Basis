function [t1,t2,opt_theta_1, opt_theta_2] = GetDegree_Total_Separate(fxy_matrix_working,gxy_matrix_working,...
    m,n,t)

global bool_preproc

% Get the degrees of polynomial f(x,y)
[r,c] = size(fxy_matrix_working);
m1 = r-1;
m2 = c-1;

% Get the degrees of polynomial g(x,y)
[r,c] = size(gxy_matrix_working);
n1 = r-1;
n2 = c-1;

% given that we know t
% it is possible that t1 and t2 take any values such that
%   1.) t1 + t2 >= t
%   2.) t1 + t2 <= 2t
%   3.) t1 <= t
%   4.) t2 <= t
%   2.) t1 <= m1
%   3.) t1 <= n1
%   4.) t2 <= m2
%   5.) t2 <= n2

%    t1 <= min(m1,n1)
%    t2 <= min(m2,n2)


minSVD = zeros( min([m1 n1 t])+1, min([m2 n2 t])+1);


lim_i = min([m1 n1 t]);
lim_j = min([m2 n2 t]);

for i = 0:1:lim_i
    count = 1;
    for j = t-i : 1 : lim_j

        
        k1 = i;
        k2 = j;
        %% Preprocessing
        switch bool_preproc
            case 'y'
                [max_mtrx_f, min_mtrx_f] = GetMaxMin(fxy_matrix_working,n1,n2,k1,k2);
                [max_mtrx_g, min_mtrx_g] = GetMaxMin(gxy_matrix_working,m1,m2,k1,k2);
                
                [theta1, theta2] = OptimalTheta(max_mtrx_f,min_mtrx_f,max_mtrx_g,min_mtrx_g);
                
            case 'n'
                theta1 = 1;
                theta2 = 1;
        end
        %% Build the Sylvester Matrix
        
        

        Cf = BuildT1(fxy_matrix_working,n1,n2,k1,k2,theta1,theta2);
        Cg = BuildT1(gxy_matrix_working,m1,m2,k1,k2,theta1,theta2);
        
        Sk1k2 = [Cf Cg];
        
        %% Get minimum Singular Value
        minSVD(i+1,j+1) = min(svd(Sk1k2));
        opt_theta_1_mtrx(i+1,j+1) = theta1;
        opt_theta_2_mtrx(i+1,j+1) = theta2;
        
        
        vec{i+j,t-j+1} = min(svd(Sk1k2));
        count = count + 1;
    end
end


vec


minSVD;