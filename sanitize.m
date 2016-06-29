function vec = sanitize(vec)

% Get number of entries in the vector
nEntries = length(vec);


% Edit - 23/02/2015 - remove any infinite values
for i = nEntries:-1:2
    if isinf(vec(i-1)) ||  vec(i-1) == 0
        vec(i-1) = vec(i);
    end
end


end
