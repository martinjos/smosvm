function check_kkt(alpha, y, K, C)

bad = zeros(4, 1);
smallval = 1e-10;
epsilon = 1e-10;

b = calcb(alpha, y, K);

for i=1:length(alpha)
    f = evalsvm(alpha, y, K, b, i);
    yfm1 = y(i) * f - 1;
    if alpha(i) < -smallval || alpha(i) > C + smallval
        % alpha out of range
        bad(1) = bad(1) + 1;
    elseif alpha(i) > smallval && alpha(i) < C - smallval
        % should have yf == 1
        t = abs(yfm1);
        if t > epsilon
            bad(3) = bad(3) + 1;
            %fprintf('Point is %g from 1\n', t);
        end
    elseif alpha(i) < C / 2
        % alpha must be close to zero - should have yf >= 1
        if yfm1 <= -epsilon
            bad(2) = bad(2) + 1;
            %fprintf('Point is %g (should be +ve)\n', yfm1);
        end
    else
        % alpha must be close to C - should have yf <= 1
        if yfm1 >= epsilon
            bad(4) = bad(4) + 1;
            %fprintf('Point is %g (should be -ve)\n', yfm1);
        end
    end
end

disp(bad);
