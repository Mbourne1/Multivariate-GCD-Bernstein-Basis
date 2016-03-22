function [x] = SolveAx_b(A,b)
%global INVERSION_METHOD


% Get the residuals associated with both QR and Pinv methods

[~,n2] = size(A);
[Q2,R] = qr(A);
R1 = R(1:n2,:);
cd = Q2'*b;
c1 = cd(1:n2,:);
x_QR = R1\c1;
res_QR = b - A*x_QR;


x_SVD = pinv(A) * b;
res_SVD = b - A*x_SVD;

if (res_SVD < res_QR)
    %winner = 'SVD';
    x = x_SVD;
else 
    %winner = 'QR';
    x = x_QR;
end


%fprintf('Minimum residual obtained by %s', winner)


end