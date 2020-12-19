function fig = stationary_sigmoid_plot(es_registered_cfpneg_logratios, ...
    es_registered_sigmoid_locations, ...
    es_registered_es_logratios, ...
    sample_names)

[~,sorted_idx] = sort(es_registered_sigmoid_locations);

colors = jet(sum(es_registered_sigmoid_locations ~= 0));
fname = strcat('../plots/stationary_sigmoid.png');
fig = figure('visible','off');
hold on;
bin_edges = -8:0.1:1;
bin_centers = bin_edges(1:length(bin_edges)-1) + 0.5 * diff(bin_edges);
legend_entries = {};
for ii=1:length(sorted_idx)
    
    idx = sorted_idx(ii);
    if isempty(es_registered_cfpneg_logratios{idx})
        continue
    end
    
    sig_loc = es_registered_sigmoid_locations(idx);
    
    % 
    counts = histcounts(es_registered_cfpneg_logratios{idx}, bin_edges);
    prob = counts ./sum(counts);
    plot(bin_centers, prob, 'Color', colors(ii,:))
    
    %
%     counts = histcounts(es_registered_es_logratios{idx}, bin_edges);
%     prob = counts ./sum(counts);
%     plot(bin_centers, prob, 'k--')
    
    plot([sig_loc, sig_loc], [0,0.05], 'Color', colors(ii,:))
    legend_entries{ii} = sample_names{idx};
end
%legend(legend_entries)
saveas(fig, fname)
end