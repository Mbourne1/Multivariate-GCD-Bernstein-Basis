function G = BuildG(t1,t2)

count =1 ;



for tot = 0:1:t1+t2+1

    for i_hat = tot:-1:0
        j_hat = tot - i_hat;
    
        
        if j_hat <= t2 && i_hat <= t1
            
            temp_vec(count) = nchoosek(t1,i_hat) * nchoosek(t2,j_hat);
            
            count = count+1;
        end
    end
end

G = diag(temp_vec);

end