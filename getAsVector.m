function f_vec = getAsVector(fxy_matrix)
% Given the polynomial f in the bernstein basis, whose coefficients are in
% matrix form, obtain the vector of the coefficients such that the order is
% increasing and the higher power of x is first.

[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

count = 1;
f_vec = zeros(r*c,1);

num_diags = r+c-1;


for tot = 0:1:num_diags
    for i = tot:-1:0
        j = tot - i;
        
        if(i > m1 || j > m2)
        else
            f_vec(count) = fxy_matrix(i+1,j+1);
            count = count + 1;
        end
        
    end
    
end


end