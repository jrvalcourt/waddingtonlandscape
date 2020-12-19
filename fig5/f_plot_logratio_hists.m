function fig = f_plot_logratio_hists(sample_name, ...
                                     logratios_3day_neg, ...
                                     logratios_3day_pos, ...
                                     logratios_es)
    fname = strcat('../plots/logratio_hists/', sample_name, '.pdf');                             
                                 
    nbins = 30;
    fig = figure('visible','off');
    hold on
    bins = linspace(-8, 3, nbins+1);
    bcenters = bins(1:nbins) + diff(bins);
    counts_neg = histcounts(logratios_3day_neg, bins);
    counts_pos = histcounts(logratios_3day_pos, bins);
    plot(bcenters, counts_neg ./ sum(counts_neg), 'k--', 'LineWidth', 2)
    plot(bcenters, counts_pos ./ sum(counts_pos), 'k', 'LineWidth', 2)
    if exist('logratios_es', 'var')
        counts_es  = histcounts(logratios_es, bins);
        plot(bcenters, counts_es ./ sum(counts_es), 'LineWidth', 2)
    end
    ylim([0,0.3])
    saveas(fig, fname)
    close()
end