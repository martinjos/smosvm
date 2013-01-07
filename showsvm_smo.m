function showsvm_smo(alpha, err, p, K, C)

smallval = 1e-10;
epsilon = 1e-100;

y = p(:,1);
points = p(:,2:3);
w = (alpha .* y)' * points;
b_old = mean(points * w' - y);

b = calcb(alpha, y, K);

fprintf('w=(%g, %g) b=%g, b_old=%g\n', w(1), w(2), b, b_old);
fprintf('dist=%g\n', b / norm(w));

axis image;
axis([0 100 0 100]);
axis manual;
hold on;

plot([0, 100], [coord(w, b, 0), coord(w, b, 100)]);

plot([0, 100], [coord(w, b + 1, 0), coord(w, b + 1, 100)], '--m');
plot([0, 100], [coord(w, b - 1, 0), coord(w, b - 1, 100)], '--m');

e = cell(length(alpha), 1);

for i=1:size(e)
    %eval = evalsvm(alpha, p(:,1), K, b, i);
    if ~kktok(i)
        e{i} = 'r';
    else
        e{i} = 'k';
    end
end

for i=1:size(p, 1)
    if p(i,1) == -1
        plot(p(i,2), p(i,3), char(['x' e{i}]));
    else
        plot(p(i,2), p(i,3), char(['o' e{i}]));
    end
    if alpha(i) > smallval && abs(i) < C - smallval
        plot(p(i,2), p(i,3), 'ks', 'MarkerSize', 20);
    end
end

function y = coord(w, b, x)

w0 = b * w / (w * w');
x = x - w0(1);
y = w0(2) - x * w(1) / w(2);
end

    function ok = kktok(i)
        %fprintf('i has length %g\n', length(i));
        %fprintf('i: %g\n', i);
        %fprintf('alpha has length %g\n', length(alpha));
        %fprintf('alpha is %g\n', alpha(i));
        %if alpha(i) < -smallval || alpha(i) > C + smallval
        %    error('alpha %d is outside 0 - C', i);
        %end
        ok = 1;
        if alpha(i) > smallval && alpha(i) < C - smallval
            if abs(err(i)) > epsilon
                ok = 0;
            end
        elseif alpha(i) <= smallval
            if y(i) * f(i) < 1 - epsilon
                ok = 0;
            end
        else
            if y(i) * f(i) > 1 + epsilon
                ok = 0;
            end
        end
    end
    function val = f(i)
        val = (alpha .* y)' * K(:,i) - b;
    end
end 