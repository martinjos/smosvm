function lotsasmosvms(down, across, p, K, y, C, epsilon)

num = down * across;

figure;
for i=1:num
    subplot(down, across, i);
    alpha = smosvm(K, y, C, epsilon);
    showsvm(alpha, p, K, true);
    obj = evalobj(alpha, y, K);
    title(obj);
end