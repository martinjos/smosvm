function obj = evalobj(alpha, y, K)

obj = sum(alpha) - (alpha .* y)' * K * (alpha .* y) / 2;