function fig = f_plot_logratio_hists_essub(sample_name, ...
                                     logratios_3day_neg, ...
                                     logratios_3day_pos, ...
                                     logratios_es)
    fname = strcat('../plots/logratio_hists_essub/', sample_name, '.png');                             
                                 
    nbins = 50;
    fig = figure('visible','off');
    hold on
    bins = linspace(-8, 2, nbins+1);
    bcenters = bins(1:nbins) + diff(bins);
    counts_neg = histcounts(logratios_3day_neg - mean(logratios_es), bins);
    counts_pos = histcounts(logratios_3day_pos - mean(logratios_es), bins);
    plot(bcenters, counts_neg ./ sum(counts_neg), 'LineWidth', 2)
    plot(bcenters, counts_pos ./ sum(counts_pos), 'LineWidth', 2)
    ylim([0,0.3])
    saveas(fig, fname)
    close()
end