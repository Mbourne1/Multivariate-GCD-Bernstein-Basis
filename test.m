function [] = test




A = sym('A', [2 2])
B = sym('B', [3 3])

%A = fliplr(triu((fliplr(A))))
B = fliplr(triu((fliplr(B))))

getInPowerBasis(A)
getInPowerBasis(B)
end