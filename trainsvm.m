function [alpha, val] = trainsvm(K, y, C, alpha)

opt = optimset('MaxFunEvals', 100000, 'MaxIter', 20000, ...
               'Algorithm', 'active-set');
%           'Display', 'off');
%               'Algorithm', 'interior-point');

if nargin < 4
    alpha = abs(randn(size(y)));
end

%disp(alpha);

[alpha, val] = fmincon(@(alpha)fun(K, y, alpha), ...
                       alpha, [], [], y', 0, zeros(size(alpha)), ...
                       C*ones(size(alpha)), [], opt);

function val = fun(K, y, alpha)

alphay = alpha .* y;

% maximise
val = - ( sum(alpha) - alphay' * K * alphay / 2 );