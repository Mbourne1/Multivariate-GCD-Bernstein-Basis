
function R1 = GetR1(Sk)
% Given the k-th Sylvester matrix get the upper triangular square matrix
% R1.
%
% Inputs.
%
% Sk

% Get QR Decomposition
% Using QR Decomposition of the sylvester matrix
[~,R] = qr(Sk);

% Take absolute values.
R = abs(R);

% Get number of rows in R1
[R1_rows,~] = size(diag(R));

% Obtain R1 the top square of the R matrix.
R1 = R(1:R1_rows,1:R1_rows);

end