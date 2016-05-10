function fw_wrt_theta2 = Differentiate_wrt_theta2(fww_matrix,th2)
% Get the partial derivative of the polynomial f(w1,w2) with respect to 
% theta_2

% Get the degree of f(w,w) with respect to y.
[~,m2] = GetDegree(fww_matrix);


temp_mat = diag((0:1:m2)./th2);

fw_wrt_theta2 =  fww_matrix * temp_mat;
    
end