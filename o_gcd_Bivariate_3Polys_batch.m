function [] = o_gcd_Bivariate_3Polys_batch()

ex_num_arr = {'1','2','3','4','5','6','7','8','9','10'};
emin_arr = {1e-8, 1e-10, 1e-12};
mean_method_arr = {'Geometric Mean Matlab Method', 'None'};
bool_alpha_theta_arr = {true,false};
low_rank_approx_method_arr = {'None','Standard STLN', 'Standard SNTLN'};
apf_method_arr = {'None'};
sylvester_build_method_arr = {'DTQ'};

parfor i1 = 1:1:length(ex_num_arr)
    ex_num = ex_num_arr{i1};
    
    for i2 = 1:1:length(emin_arr)
        
        emin = emin_arr{1};
        emax = 1e-12;
        
        for mean_method_elem = mean_method_arr
            
            mean_method = mean_method_elem{1};
            
            for bool_alpha_theta_elem = bool_alpha_theta_arr
                bool_alpha_theta = bool_alpha_theta_elem{1};
                
                for low_rank_approx_method_elem = low_rank_approx_method_arr
                    low_rank_approx_method = low_rank_approx_method_elem{1};
                    
                    for apf_method_elem = apf_method_arr
                        apf_method =  apf_method_elem{1};
                        
                        for sylvester_build_method_elem = sylvester_build_method_arr
                            sylvester_matrix_variant = sylvester_build_method_elem{1};
                            
                            close all; clc;
                            
                            my_filename = 'log_GCD_3Polys.txt'
                            
                            try
                                
                                o_gcd_Bivariate_3Polys(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, sylvester_matrix_variant);
                                fileId = fopen(my_filename, 'a')
                                fprintf(fileId,'%s %s \n',datetime('now'), 'Success');
                                fclose(fileId);
                                
                            catch err
                                
                                fileId = fopen(my_filename, 'a')
                                fprintf(fileId,'%s %s \n\n\n\',datetime('now'), getReport(err));
                                fclose(fileId);
                                
                            end
                        end
                    end
                end
                
            end
        end
    end
end


end