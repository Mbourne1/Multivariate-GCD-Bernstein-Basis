function [] = plotR1Diagonals_1Dimensional(arr_R1, limits_k1, limits_t1)
%
% % Inputs
%
% arr_R1 : (Array)
%
% limits_k1 : ([Int Int])
%
% limits_t1 : (Int Int)

% Get upper and lower limits
myLowerLimit_t1 = limits_k1(1);
myUpperLimit_t1 = limits_k1(2);

lowerLimit_t1 = limits_t1(1);
upperLimit_t1 = limits_t1(2);

% Get number of Subresultants
nSubresultants_t1 = myUpperLimit_t1 - myLowerLimit_t1 +1 ;

global SETTINGS

figure_name = sprintf('Diagonals of R1 from QR decomp of %s',SETTINGS.SYLVESTER_MATRIX_VARIANT);
figure('name', figure_name)
hold on

for i1 = 1:1:nSubresultants_t1
  
        k1 = myLowerLimit_t1 + (i1 -1);
  
        temp_vec = diag(arr_R1{i1, i2});
        
        vec_k1 = k1.*ones(length(temp_vec), 1);
              
        plot(vec_k1, temp_vec,'*');
        
end

xlabel('k_{1}')

ylabel('log_{10} Row Diagonals of R_{1}')
vline(lowerLimit_t1);
vline(upperLimit_t1);
grid on
hold off

end