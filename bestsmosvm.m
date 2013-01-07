function bestalpha = bestsmosvm(K, y, C, epsilon, num)

bestobj = -Inf;
bestalpha = [];

for i=1:num
    alpha = smosvm(K, y, C, epsilon);
    obj = evalobj(alpha, y, K);
    if obj > bestobj
        bestobj = obj;
        bestalpha = alpha;
    end
    fprintf('Done %d\n', i);
end