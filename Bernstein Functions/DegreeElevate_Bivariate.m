function deg_elv_fxy_matrix = DegreeElevate_Bivariate(fxy, ad_rows, ad_cols)
% Degree elevate the polynomial f(x,y) so that f(x,y) is a polynomial with
% basis elements
% B_{m1 + ad_rows}(x) in x
% B_{m2 + ad_cols}(y) in y
%
%
% % Inputs:
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% ad_rows : (Int) Number of degree elevations with respect to x
%
% ad_cols : (Int) Number of degree elevations with respect to y.



% let r_f be the number of rows being added
r_f = ad_rows;

% let s_f be the number of cols being added
s_f = ad_cols;

% Degree elevate fxy_matrix
%deg_elv_fxy_matrix = zeros(m1 + ad_rows, m2 + ad_cols);

gxy = ones(r_f + 1, s_f + 1);
deg_elv_fxy_matrix = Bernstein_Multiply_Bivariate(fxy,gxy);


% % for each column j
% for j = 0:1:m2 + s_f
% 
%     % for each row i
%     for i = 0:1:m1 + r_f
% 
%         deg_elv_fxy_matrix(i+1,j+1) = summing_function(fxy, m1, m2, r_f, s_f, i, j);
%     end
% 
% end

%deg_elv_fxy_matrix_test - deg_elv_fxy_matrix

end

% function sum = summing_function(fxy, m1, m2, r, s, i, j)
% %
% % % Inputs
% %
% % fxy : (Matrix) Coefficients of polynomial f(x,y)
% %
% % m1 : (Int) Degree of polynomial f(x,y) with respect to x
% %
% % m1 : (Int) Degree of polynomial f(x,y) with respect to y
% %
% % r : (Int) Number of degree elevations in x
% %
% % s : (Int) Number of degree elevations in y
% %
% % i :
% %
% % j :
% 
% % Initialise a sum
% sum = 0;
% 
% for k = max(0,i-r) : 1 : min(m1,i)
%     for l = max(0,j-s) : 1 : min(m2,j)
%         sum = sum +...
%             (...
%             fxy(k+1,l+1) *...
%             nchoosek(m1,k) *...
%             nchoosek(r,i-k) *...
%             nchoosek(m2,l) *...
%             nchoosek(s,j-l) /...
%             ...
%             (nchoosek(m1+r,i) *...
%             nchoosek(m2+s,j))...
%             );
%     end
% end
% 
% end

