function value = evalsvm(alpha, y, K, b, i)

value = K(i,:) * (alpha .* y) - b;