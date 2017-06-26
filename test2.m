
for m1 = 1:1:7

prod = 1;

for i1 = 0:1:m1
   
    prod = prod .* nchoosek(m1,i1);
    
end

display(prod)

vec(i1+1) = prod;

end

display(vec)