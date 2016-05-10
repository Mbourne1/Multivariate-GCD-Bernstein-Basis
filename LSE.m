function y = LSE(~,f,C,g)
% This function uses the QR decomposition to solve the LSE problem
% minimise ||Ey-f||  subject to  Cy=g
% The output is the vector y.


    [m1,~] = size(C);
    
    [Q,R] = qr(C'); 
    
    R1 = R(1:m1,:);
    
    w1 = R1'\g;
    
    %Q1 = Q(:,1:m1);
    Q2 = Q(:,m1+1:end);

    w2 = Q2'*f;

    y = Q*[w1;w2];



end

