function [arr_hxy] = Deconvolve_Bivariate_Batch(arr_fxy, vDeg_fxy)
% Get the set of polynomials h_{i} given by the deconvolution of the
% polynomials f_{i}, where h_{i} = f_{i-1}/f_{i}
%
% % Inputs.
%
%
% arr_fxy : (Array of Matrices) Array of polynomials f_{i}(x,y)
%
% vDeg_fxy : (Vector) Vector containing the total degree of the polynomials f{i}
%
%
% % Outputs.
%
%
% arr_hxy : (Array of Matrices) Array of polynomials h{i}(x,y)

% %
% %
% Form the left hand side matrix

% Get number of polynomials in the array arr_fxy of f_{i}(x,y)
nPolys_arr_fxy = size(arr_fxy, 1);

% Get the degrees of polynomials f_{i}(x,y)
vDeg_x_fxy = zeros(nPolys_arr_fxy, 1);
vDeg_y_fxy = zeros(nPolys_arr_fxy, 1);
for i = 1:1:nPolys_arr_fxy
    
    [vDeg_x_fxy(i), vDeg_y_fxy(i)]  = GetDegree_Bivariate(arr_fxy{i});
     
end


% % Preprocess polynomials f_{i,j}(x,y)
global SETTINGS
if( SETTINGS.PREPROC_DECONVOLUTIONS)
    
    [th1, th2] = GetOptimalTheta(arr_fxy);
   
    
else
    
    th1 = 1;
    th2 = 1;
    
end


arr_fww = GetPolynomialArrayWithThetas(arr_fxy, th1, th2);

% Build the matrix C(f1,...,fd)
C_fw = BuildC(arr_fww);


% Form the right hand side vector
vRHS = BuildRHS_vec(arr_fww);


% Get vector of coefficients of the polynomials h_{i}(x,y)
v_hww = SolveAx_b(C_fw, vRHS);

% split solution vector into polynomials h_{i}(x,y)
%vDeg_hxy = vDeg_fxy(1 : end-1) - vDeg_fxy(2 : end);

vDeg_x_hxy = vDeg_x_fxy(1 : end - 1) - vDeg_x_fxy(2 : end);
vDeg_y_hxy = vDeg_y_fxy(1 : end - 1) - vDeg_y_fxy(2 : end);

% Get array of polynomials
arr_hww = GetPolynomialArrayFromVector(v_hww, vDeg_x_hxy, vDeg_y_hxy);

% Get without thetas
arr_hxy = GetPolynomialArrayWithoutThetas(arr_hww, th1, th2);


end

function C_fxy = BuildC(arr_fxy)
% Build the matrix C(f1,...,fd)
%
% Inputs.
%
% arr_fxy : (Array of Matrices) Array of polynomials f(x,y)
%
% Outputs.
%
% C_fxy : (Matrxi)

% Get number of polynomials in f_{i}(x,y)
nPolys_arr_fxy = length(arr_fxy);

% Get degree of each polynomial f_{i}(x,y)
vDeg_x_fxy = zeros(nPolys_arr_fxy, 1);
vDeg_y_fxy = zeros(nPolys_arr_fxy, 1);

for i = 1 : 1 : nPolys_arr_fxy
    
    [vDeg_x_fxy(i), vDeg_y_fxy(i)] = GetDegree_Bivariate(arr_fxy{i});
     
end

arr_DT1Q1 = cell(nPolys_arr_fxy, 1);

% For each of the polynomials excluding the first f_{1},...,f_{d}
for i = 2 : 1 : nPolys_arr_fxy
    
    % Get the degree of f{i-1}
    n1 = vDeg_x_fxy(i - 1);
    n2 = vDeg_y_fxy(i - 1);
    
    % Get the degree of f{i}
    m1 = vDeg_x_fxy(i);
    m2 = vDeg_y_fxy(i);
    
    % Temporarily call the ith entry f(x,y)
    fxy = arr_fxy{i};
    
    % Build the matrix T_{n-m}(f(x,y))
    
    D = BuildD_Bivariate_2Polys(m1, m2, n1 - m1, n2 - m2);
    T1 = BuildT1_Bivariate(fxy, n1 - m1, n2 - m2);
    Q1 = BuildQ1_Bivariate(n1 - m1, n2 - m2);
   
    arr_DT1Q1{i-1} = D*T1*Q1;
    
end



C_fxy = blkdiag(arr_DT1Q1{:});


end


