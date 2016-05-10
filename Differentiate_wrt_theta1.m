function fw_wrt_theta1 = Differentiate_wrt_theta1(fww_matrix,th1)
% Get the partial derivative of the polynomial f(w,w) with repsect to
% theta_{1}

% Get the degree of f(w,w) with respect to x
[m1,~] = GetDegree(fww_matrix);


temp_mat = diag((0:1:m1)./th1);

fw_wrt_theta1 = temp_mat * fww_matrix;

end