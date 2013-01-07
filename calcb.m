function b = calcb(alpha, y, K)

%svs = find(alpha > 1e-8);
svs = find(alpha ~= 0); % works with smosvm
b = mean(sum(repmat(alpha .* y, 1, length(svs)) .* K(:,svs)) - y(svs)');