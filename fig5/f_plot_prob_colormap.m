function fig = f_plot_prob_colormap(sample_name, sigmoidA, sigmoidB)
lims = [-8,3];

fig = figure('visible','off');
nbins = 10000;
fname = strcat('../plots/prob_colormap/', sample_name, '.png');
bins = linspace(lims(1), lims(2), nbins+1);
bcenters = bins(1:nbins) + 0.5 .* diff(bins);
image(lims, [0,1], f_sigmoid(bcenters, sigmoidA, sigmoidB),...
      'CDataMapping','scaled')

temp = yellow_blue_cmap();
colormap(temp);
saveas(fig, fname)
close()
end