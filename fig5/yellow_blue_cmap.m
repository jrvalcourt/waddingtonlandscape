function cmap = yellow_blue_cmap()
ncolors = 1000;
x = 1:ncolors;
x = x ./ ncolors;

% yellow 209 171 43
% blue 40 126 194
blues_r =  40./255 + x .* (255-40)  ./ 255;
blues_g = 126./255 + x .* (255-126) ./ 255;
blues_b = 194./255 + x .* (255-194) ./ 255;
yellow_r = 1 - x .* (255-209) ./ 255;
yellow_g = 1 - x .* (255-171) ./ 255;
yellow_b = 1 - x .* (255-43)  ./ 255;
blues   = [blues_r;  blues_g;  blues_b]';
yellows = [yellow_r; yellow_g; yellow_b]';
cmap = [blues; yellows];

end

