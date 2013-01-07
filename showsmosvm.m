function [alpha, b] = showsmosvm(p, K, y, C, epsilon)
    figure;
    function plotfunc(code, arg)
        if isequal(code, 'do')
            plot(p(arg(1),2), p(arg(1),3), 'or', 'MarkerSize', 30);
            plot(p(arg(2),2), p(arg(2),3), 'ob', 'MarkerSize', 30);
        elseif isequal(code, 'done')
            clf;
            showsvm_smo(arg{1}, arg{2}, p, K, C);
        end
        pause(1);
    end
    [alpha, b] = smosvm_gfx(K, y, C, epsilon, @plotfunc);
end