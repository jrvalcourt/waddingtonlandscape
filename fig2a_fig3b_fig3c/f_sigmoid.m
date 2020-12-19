function y = f_sigmoid(x,a,b)
    y = 1 ./ (1 + exp(-a .* (x - b)));
end