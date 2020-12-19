function fig = f_plot_sigmoid_shift(sample_name, sigmoidA, ...
    sigmoidB_neg, sigmoidB_pos);


lims = [-8,3];

fig = figure();
nbins = 1000;
fname = strcat('../plots/sigmoid_shift/', sample_name, '.pdf');
bins = linspace(lims(1), lims(2), nbins+1);
bcenters = bins(1:nbins) + 0.5 .* diff(bins);
hold off
if length(sigmoidB_pos) > 1
    mean_sigmoidB_pos = mean(sigmoidB_pos);
    std_sigmoidB_pos  = std(sigmoidB_pos);
    patch([bcenters fliplr(bcenters)], ...
        [f_sigmoid(bcenters, sigmoidA, ...
                   mean_sigmoidB_pos - std_sigmoidB_pos) ...
        fliplr(f_sigmoid(bcenters, sigmoidA, ...
                         mean_sigmoidB_pos + std_sigmoidB_pos))], ...
        [0.9, 0.9, 0.9], ...
        'EdgeColor', 'none')
    hold on
    plot(bcenters, f_sigmoid(bcenters, sigmoidA, mean_sigmoidB_pos), ...
        'k', 'LineWidth', 2)
else
    plot(bcenters, f_sigmoid(bcenters, sigmoidA, sigmoidB_pos(1)), ...
         'k', 'LineWidth', 2)
    hold on
end
plot(bcenters, f_sigmoid(bcenters, sigmoidA, sigmoidB_neg), 'k--', ...
    'LineWidth', 2)
% for ii=1:length(sigmoidB_pos)
%     plot(bcenters, f_sigmoid(bcenters, sigmoidA, sigmoidB_pos(ii)), ...
%         'k', 'LineWidth', 2)
% end
saveas(fig, fname)
close()
end