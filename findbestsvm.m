function save = findbestsvm(K, y, C, num)

save = zeros(num, 1 + length(y));

for i=1:num
    [alpha, val] = trainsvm(K, y, C);
    save(i,:) = [val, alpha'];
end