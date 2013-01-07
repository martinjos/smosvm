function table = testsvm(alpha, y, K, b)

table = zeros(2, 2);
for i=1:length(y)
    newy = sign(evalsvm(alpha, y, K, b, i));
    same_idx = 1;
    if newy ~= y(i)
        same_idx = 2;
    end
    class_idx = 1;
    if y(i) == 1
        class_idx = 2;
    end
    table(class_idx, same_idx) = table(class_idx, same_idx) + 1;
end