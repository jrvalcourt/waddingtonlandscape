function [figs, results] = f_results_from_mat_jim_version(gene, label, dataDir_3day, dataDir_endpoint, baseName_3day, baseName_endpoint, numBins)
% 3 day analysis: KL

figs = [];
results = {};
results{1} = gene;
results{2} = label;

% Load 3-day data
Data = load([dataDir_3day baseName_3day]);
Data = Data.Data;
% Data is a cell array of the following order:
    % Cell 1: CFP_negative OCT4
    % Cell 2: CFP_negative SOX2
    % Cell 3: CFP_positive OCT4
    % Cell 4: CFP_positive SOX2
% Turn this cell array into two tables
varNames = {'Red', 'Green'};
DATA_neg_full = table(Data{1}, Data{2}, 'VariableNames', varNames);
DATA_pos_full = table(Data{3}, Data{4}, 'VariableNames', varNames);

% Find cells with values above 0 in both columns (CFP-negative)
red_above_idx = find(DATA_neg_full.Red > 0);
green_above_idx = find(DATA_neg_full.Green > 0);
all_above_idx = intersect(red_above_idx, green_above_idx);
DATA_neg = DATA_neg_full(all_above_idx,:);
DATA_neg_3day = DATA_neg_full(all_above_idx,:);

% Find cells with positive values in both columns (CFP-positive)
red_above_idx = find(DATA_pos_full.Red > 0);
green_above_idx = find(DATA_pos_full.Green > 0);
all_above_idx = intersect(red_above_idx, green_above_idx);
DATA_pos = DATA_pos_full(all_above_idx,:);
DATA_pos_3day = DATA_pos_full(all_above_idx,:);

% Countour plots: 3day
n_grid = 30;
sox2_min = 10;
sox2_max = 2 .* 10^5;
oct4_min = 10;
oct4_max = 10^5;

figs(1) = figure();
centers = {10.^linspace(log10(oct4_min),log10(oct4_max), n_grid), ...
           10.^linspace(log10(sox2_min),log10(sox2_max), n_grid)};
[n,c] = hist3([DATA_neg.Red, DATA_neg.Green], centers);
contour(c{1}, c{2}, n', 10, 'LineWidth', 2)
ylabel('SOX2::YFP (a.u.)','FontSize',10);
xlabel('OCT4::YFP (a.u.)','FontSize',10);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([oct4_min, oct4_max])
ylim([sox2_min, sox2_max])
hold off

figs(2) = figure();
centers = {10.^linspace(log10(oct4_min),log10(oct4_max), n_grid), ...
           10.^linspace(log10(sox2_min),log10(sox2_max), n_grid)};
[n,c] = hist3([DATA_pos.Red, DATA_pos.Green], centers);
contour(c{1},c{2},n',10, 'LineWidth', 2)
ylabel('SOX2::YFP (a.u.)','FontSize',10);
xlabel('OCT4::YFP (a.u.)','FontSize',10);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([oct4_min, oct4_max])
ylim([sox2_min, sox2_max])
hold off

% Take the log of each data point
neg_log_red   = log(DATA_neg.Red);
neg_log_green = log(DATA_neg.Green);
pos_log_red   = log(DATA_pos.Red);
pos_log_green = log(DATA_pos.Green);

% Take the ratio
DATA_neg.LogRatio = log2(DATA_neg.Red ./ DATA_neg.Green);
DATA_pos.LogRatio = log2(DATA_pos.Red ./ DATA_pos.Green);

% Endpoint analysis: 

Data = load([dataDir_endpoint baseName_endpoint]);
Data = Data.Data;
    % Cell 1: CFP_negative OCT4
    % Cell 2: CFP_negative SOX2
    % Cell 3: CFP_positive OCT4
    % Cell 4: CFP_positive SOX2
varNames = {'Red','Green'};
DATA_neg_full = table(Data{1}, Data{2}, 'VariableNames', varNames);
DATA_pos_full = table(Data{3}, Data{4}, 'VariableNames', varNames);

% Find cells with values above 0 in both columns (CFP-negative)
red_above_idx = find(DATA_neg_full.Red > 0);
green_above_idx = find(DATA_neg_full.Green > 0);
all_above_idx = intersect(red_above_idx, green_above_idx);
DATA_neg = DATA_neg_full(all_above_idx,:);

% Find cells with positive values in both columns (CFP-positive)
red_above_idx = find(DATA_pos_full.Red > 0);
green_above_idx = find(DATA_pos_full.Green > 0);
all_above_idx = intersect(red_above_idx, green_above_idx);
DATA_pos = DATA_pos_full(all_above_idx,:);

% Endpoint figure
n_grid = 30;
sox2_min = 50;
sox2_max = 2 * 10^5;
oct4_min = 10;
oct4_max = 5 * 10^4;

figs(3) = figure();
centers = {10.^linspace(log10(oct4_min), log10(oct4_max), n_grid), ...
           10.^linspace(log10(sox2_min), log10(sox2_max), n_grid)};
[n,c] = hist3([DATA_neg.Red, DATA_neg.Green], centers);
contour(c{1}, c{2}, n', 10, 'LineWidth', 2)
ylabel('SOX2::YFP (a.u.)','FontSize',10);
xlabel('OCT4::YFP (a.u.)','FontSize',10);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([oct4_min, oct4_max])
ylim([sox2_min, sox2_max])
hold off

figs(4) = figure();
centers = {10.^linspace(log10(oct4_min), log10(oct4_max), n_grid), ...
           10.^linspace(log10(sox2_min), log10(sox2_max), n_grid)};
[n,c] = hist3([DATA_pos.Red, DATA_pos.Green], centers);
contour(c{1}, c{2}, n', 10, 'LineWidth', 2)
ylabel('SOX2::YFP (a.u.)','FontSize',10);
xlabel('OCT4::YFP (a.u.)','FontSize',10);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([oct4_min, oct4_max])
ylim([sox2_min, sox2_max])
hold off

% Gating 

poly_file_be = ['polys/', label, '_BE.mat'];
poly_file_me = ['polys/', label, '_ME.mat'];

if exist(poly_file_be, 'file')
    load(poly_file_be, 'coords_BE')
else
    hold off
%     n_grid = 30;
%     sox2_min = 10^2;
%     sox2_max = 1 * 10^5;
%     oct4_min = 10;
%     oct4_max = 2 * 10^5;
    % Visualize a manageable subset of CFP-negative cells
    if height(DATA_neg) > 50000
        rand_idx = randsample(height(DATA_neg), 50000);
        DATA_neg_disp = DATA_neg(rand_idx,:);
    else
        DATA_neg_disp = DATA_neg;
    end
    % Plot and gate endpoints of CFP-negative populations
    scatter(DATA_neg_disp.Red, DATA_neg_disp.Green, '.', 'MarkerEdgeColor','black', 'MarkerFaceAlpha', 0.03,'MarkerEdgeAlpha', 0.03);
    hold on
    scatter(DATA_pos.Red, DATA_pos.Green, '.', 'MarkerEdgeColor',[0 0.5 1], 'MarkerFaceAlpha', 0.03,'MarkerEdgeAlpha', 0.03);
    centers = {10.^linspace(log10(oct4_min),log10(oct4_max),n_grid), 10.^linspace(log10(sox2_min),log10(sox2_max),n_grid)};
    [n,c] = hist3([DATA_neg.Red, DATA_neg.Green], centers);
    contour(c{1},c{2},n',10, 'LineColor', 'black', 'LineWidth', 0.5)
    [n2,c2] = hist3([DATA_pos.Red, DATA_pos.Green], centers);
    contour(c2{1},c2{2},n2',10, 'LineStyle', '--', 'LineColor', [0 0.5 1], 'LineWidth', 1)
    title('Draw bipotent ectoderm polygon:')
    ylabel('SOX2::YFP','FontSize',10);
    xlabel('OCT4::YFP','FontSize',10);
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    h = drawpolygon();
    coords_BE = h.Position;
    save(poly_file_be, 'coords_BE')
end
hold off 

if exist(poly_file_me, 'file')
    load(poly_file_me, 'coords_ME')
else
    hold off
%     n_grid = 30;
%     sox2_min = 10^2;
%     sox2_max = 1 * 10^5;
%     oct4_min = 10;
%     oct4_max = 2 * 10^5;
    % Visualize a manageable subset of CFP-negative cells
    if height(DATA_neg) > 50000
        rand_idx = randsample(height(DATA_neg), 50000);
        DATA_neg_disp = DATA_neg(rand_idx,:);
    else
        DATA_neg_disp = DATA_neg;
    end
    % Plot and gate endpoints of CFP-negative populations
    scatter(DATA_neg_disp.Red, DATA_neg_disp.Green, '.', ...
        'MarkerEdgeColor','black', 'MarkerFaceAlpha', 0.03,...
        'MarkerEdgeAlpha', 0.03);
    hold on
    scatter(DATA_pos.Red, DATA_pos.Green, '.', ...
        'MarkerEdgeColor',[0 0.5 1], 'MarkerFaceAlpha', 0.03,...
        'MarkerEdgeAlpha', 0.03);
    centers = {10.^linspace(log10(oct4_min),log10(oct4_max),n_grid), ...
               10.^linspace(log10(sox2_min),log10(sox2_max),n_grid)};
    [n,c] = hist3([DATA_neg.Red, DATA_neg.Green], centers);
    contour(c{1},c{2},n',10, 'LineColor', 'black', 'LineWidth', 0.5)
    [n2,c2] = hist3([DATA_pos.Red, DATA_pos.Green], centers);
    contour(c2{1},c2{2},n2',10, 'LineStyle', '--', ...
        'LineColor', [0 0.5 1], 'LineWidth', 1)
    title('Draw mesendoderm polygon:')
    ylabel('SOX2::YFP','FontSize',10);
    xlabel('OCT4::YFP','FontSize',10);
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    h = drawpolygon();
    coords_ME = h.Position;
    save(poly_file_me, 'coords_ME')
end
hold off 

% Index BE cells in both CFP+ and CFP- populations
cfp_neg_BE = inpolygon(DATA_neg.Red, DATA_neg.Green, ...
                       coords_BE(:,1), coords_BE(:,2));
cfp_neg_ME = inpolygon(DATA_neg.Red, DATA_neg.Green, ...
                       coords_ME(:,1), coords_ME(:,2));                   
cfp_pos_BE = inpolygon(DATA_pos.Red, DATA_pos.Green, ...
                       coords_BE(:,1), coords_BE(:,2));
cfp_pos_ME = inpolygon(DATA_pos.Red, DATA_pos.Green, ...
                       coords_ME(:,1), coords_ME(:,2)); 
                   
ratios_pos = DATA_pos_3day.Red ./ DATA_pos_3day.Green;
ratios_neg = DATA_neg_3day.Red ./ DATA_neg_3day.Green;

sorted_ratios_pos = sort(ratios_pos);
sorted_ratios_neg = sort(ratios_neg);

frac_BE_pos = sum(cfp_pos_BE) ./ (sum(cfp_pos_BE) + sum(cfp_pos_ME));
frac_BE_neg = sum(cfp_neg_BE) ./ (sum(cfp_neg_BE) + sum(cfp_neg_ME));

ratio_pos = sorted_ratios_pos(floor(frac_BE_pos * length(sorted_ratios_pos)));
ratio_neg = sorted_ratios_neg(floor(frac_BE_neg * length(sorted_ratios_neg)));

% Countour plots: 3days with lines
n_grid = 30;
sox2_min = 10;
sox2_max = 2 .* 10^5;
oct4_min = 10;
oct4_max = 10^5;

figs(5) = figure();
centers = {10.^linspace(log10(oct4_min), log10(oct4_max), n_grid), ...
           10.^linspace(log10(sox2_min), log10(sox2_max), n_grid)};
[n,c] = hist3([DATA_neg_3day.Red, DATA_neg_3day.Green], centers);
contour(c{1}, c{2}, n', 10, 'LineWidth', 2)
hold on
plot([oct4_min,              oct4_max], ...
     [oct4_min ./ ratio_neg, oct4_max ./ ratio_neg], 'k--')
plot([oct4_min,              oct4_max], ...
     [oct4_min ./ ratio_pos, oct4_max ./ ratio_pos], 'b--')
ylabel('SOX2:YFP (a.u.)','FontSize',10);
xlabel('OCT4:YFP (a.u.)','FontSize',10);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([oct4_min, oct4_max])
ylim([sox2_min, sox2_max])
hold off

figs(6) = figure();
centers = {10.^linspace(log10(oct4_min), log10(oct4_max), n_grid), ...
           10.^linspace(log10(sox2_min), log10(sox2_max), n_grid)};
[n,c] = hist3([DATA_pos_3day.Red, DATA_pos_3day.Green], centers);
contour(c{1}, c{2}, n', 10, 'LineWidth', 2)
hold on
plot([oct4_min,              oct4_max], ...
     [oct4_min ./ ratio_neg, oct4_max ./ ratio_neg], 'k--')
plot([oct4_min,              oct4_max], ...
     [oct4_min ./ ratio_pos, oct4_max ./ ratio_pos], 'b--')
ylabel('SOX2:YFP (a.u.)','FontSize',10);
xlabel('OCT4:YFP (a.u.)','FontSize',10);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([oct4_min, oct4_max])
ylim([sox2_min, sox2_max])
hold off

figs(7) = figure();
bins = -5:0.1:1;
[n_neg,bins_neg] = hist(log10(ratios_neg), bins);
[n_pos,bins_pos] = hist(log10(ratios_pos), bins);
plot(bins_neg, n_neg ./ sum(n_neg), 'ko-')
hold on
text_position = bins_neg(3);
text(text_position, 0.25, strcat("frac BE neg = ", num2str(frac_BE_neg)))
text(text_position, 0.2, strcat("frac BE pos = ", num2str(frac_BE_pos)))
plot(bins_pos, n_pos ./ sum(n_pos), 'bo-')
plot([log10(ratio_neg), log10(ratio_neg)], [0,0.3], 'k--')
plot([log10(ratio_pos), log10(ratio_pos)], [0,0.3], 'b--')
hold off

figs(8) = figure();
bins = -5:0.1:1;
[n_neg,bins_neg] = hist(log10(ratios_neg), bins);
[n_pos,bins_pos] = hist(log10(ratios_pos), bins);
neg_cumsum = cumsum(n_neg);
neg_cumsum = neg_cumsum ./ max(neg_cumsum);
plot(bins_neg, neg_cumsum, 'ko-')
hold on
text_position = bins_neg(3);
text(text_position, 0.9, strcat("frac BE neg = ", num2str(frac_BE_neg)))
text(text_position, 0.85, strcat("frac BE pos = ", num2str(frac_BE_pos)))
pos_cumsum = cumsum(n_pos);
pos_cumsum = pos_cumsum ./ max(pos_cumsum);
plot(bins_pos, pos_cumsum, 'bo-')
plot([log10(ratio_neg), log10(ratio_neg)], [0,1], 'k--')
plot([log10(ratio_pos), log10(ratio_pos)], [0,1], 'b--')
delta_cutoff = log10(ratio_neg) - log10(ratio_pos);
diff_idx = find(neg_cumsum >= frac_BE_neg, 1, 'first');
diffs = pos_cumsum - neg_cumsum;
actual_diff  = frac_BE_pos - frac_BE_neg;
realistic_diff = pos_cumsum(diff_idx) - neg_cumsum(diff_idx);
if realistic_diff > 0
    biggest_diff = max(diffs);
else
    biggest_diff = min(diffs);
end
hold off

figs(9) = figure();
bins = 0:0.01:1;
[n_neg,bins_neg] = hist(ratios_neg, bins);
[n_pos,bins_pos] = hist(ratios_pos, bins);
neg_cumsum = cumsum(n_neg);
neg_cumsum = neg_cumsum ./ max(neg_cumsum);
plot(bins_neg, neg_cumsum, 'ko-')
hold on
text_position = bins_neg(3);
text(text_position, 0.9, strcat("frac BE neg = ", num2str(frac_BE_neg)))
text(text_position, 0.85, strcat("frac BE pos = ", num2str(frac_BE_pos)))
pos_cumsum = cumsum(n_pos);
pos_cumsum = pos_cumsum ./ max(pos_cumsum);
plot(bins_pos, pos_cumsum, 'bo-')
plot([ratio_neg, ratio_neg], [0,1], 'k--')
plot([ratio_pos, ratio_pos], [0,1], 'b--')
delta_cutoff = ratio_neg - ratio_pos;
diff_idx = find(neg_cumsum >= frac_BE_neg, 1, 'first');
diffs = pos_cumsum - neg_cumsum;
actual_diff  = frac_BE_pos - frac_BE_neg;
realistic_diff = pos_cumsum(diff_idx) - neg_cumsum(diff_idx);
if realistic_diff > 0
    biggest_diff = max(diffs);
else
    biggest_diff = min(diffs);
end
hold off

% Calculate BE proportions
BE_neg = sum(cfp_neg_BE) / (sum(cfp_neg_BE) + sum(cfp_neg_ME));
BE_pos = sum(cfp_pos_BE) / (sum(cfp_pos_BE) + sum(cfp_pos_ME));
%[gene, ' ', label, ' ', num2str(1-BE_neg), ' ', num2str(1-BE_pos)]

% Create histogram data

% Take the ratio of the fluorescence levels
DATA_neg_3day.Ratio = DATA_neg_3day.Red ./ DATA_neg_3day.Green;
DATA_pos_3day.Ratio = DATA_pos_3day.Red ./ DATA_pos_3day.Green;

% Take the log of each data point and add a pseudocount to deal with zeros
pseudocount = 0.00001;
DATA_neg_3day.LogRatio = log(DATA_neg_3day.Ratio + pseudocount);
DATA_pos_3day.LogRatio = log(DATA_pos_3day.Ratio + pseudocount);
h = histogram(DATA_neg_3day.LogRatio, numBins);
hist_neg = h.Values;
negEdges = h.BinEdges;
h2 = histogram(DATA_pos_3day.LogRatio, 'BinEdges', negEdges);
hist_pos = h2.Values;

% calc KL
results{3} = f_kldiv(hist_pos, hist_neg, 0.000001);

% Calculate log of ratios
results{4} = log2( (BE_pos / (1 - BE_pos)) ...
                 / (BE_neg / (1 - BE_neg)) );
             
results{5} = delta_cutoff;
results{6} = [biggest_diff, actual_diff, realistic_diff];

% find the sigmoid shift that would give us this fraction mesendoderm
sigmoidA = 8.37755655; % TODO make this pull the most up-to-date value 
                       % from the pickle file
results{7} = [f_infer_sigmoid_loc(log(ratios_neg)/log(2), 1 - BE_neg, sigmoidA), ...
              f_infer_sigmoid_loc(log(ratios_pos)/log(2), 1 - BE_pos, sigmoidA)];

end

function kl = f_kldiv(query, ref, pseudocount)
    
    % first make sure that the two input histograms have the same
    % dimensions since the probability distributions MUST have support over
    % the same domain
    dims_query = size(query);
    dims_ref   = size(ref);
    assert(all(dims_query == dims_ref));
    
    % get a total count to normalize each 2D histogram to estimate prob
    total_count_query = 0;
    total_count_ref   = 0;
    for ii=1:dims_query(1)
        for jj=1:dims_query(2)
            total_count_query = total_count_query + query(ii,jj) + pseudocount;
            total_count_ref   = total_count_ref   + ref(ii,jj)   + pseudocount;
        end
    end
    
    % calculate the kl divergence
    kl = 0;
    for ii=1:dims_query(1)
        for jj=1:dims_query(2)
            p = (query(ii,jj) + pseudocount) / total_count_query;
            q = (ref(ii,jj)   + pseudocount) / total_count_ref;
            kl = kl + (p * log2(p / q));
        end
    end
end