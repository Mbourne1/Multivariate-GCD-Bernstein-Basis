function [] = plotR1Diagonals(arr_R1, my_limits_t1, my_limits_t2)
%
% % Inputs
%
% arr_R1 : (Array)
%
% my_limits_t1 : ([Int Int])
%
% my_limits_t2 : ([Int Int])

% Get upper and lower limits
lowerLimit_t1 = my_limits_t1(1);
upperLimit_t1 = my_limits_t1(2);
lowerLimit_t2 = my_limits_t2(1);
upperLimit_t2 = my_limits_t2(2);

% Get number of Subresultants
nSubresultants_t1 = upperLimit_t1 - lowerLimit_t1 +1 ;
nSubresultants_t2 = upperLimit_t2 - lowerLimit_t2 +1 ;

global SETTINGS

figure_name = sprintf('Diagonals of R1 from QR decomp of %s',SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name', figure_name)
hold on

for i1 = 1:1:nSubresultants_t1
    for i2 = 1:1:nSubresultants_t2
        
        k1 = lowerLimit_t1 + (i1 -1);
        k2 = lowerLimit_t2 + (i2 -1);
        
        temp_vec = diag(arr_R1{i1, i2});
        
        v_k1 = k1.*ones(length(temp_vec), 1);
        v_k2 = k2.*ones(length(temp_vec), 1);
        
        plot3(v_k1 ,v_k2, log10(temp_vec), '*');
        
    end
end
xlabel('k_{1}')
ylabel('k_{2}')
zlabel('log_{10} Row Diagonals of R_{1}')

grid on
hold off

end