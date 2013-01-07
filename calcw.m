function w = calcw(alpha, y, X)

disp(size(X));
w = (alpha .* y)' * X;