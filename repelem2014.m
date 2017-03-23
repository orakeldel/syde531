function res = repelem2014(a, b) % repeats elements

index = zeros(1, sum(b));
index([1 cumsum(b(1:end-1))+1]) = 1;
res = a(cumsum(index));

end