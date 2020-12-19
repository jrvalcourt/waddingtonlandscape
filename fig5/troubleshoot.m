clear all
close all

refList = readtable('reference_list.xlsx');
dataDir_3day = '../../data/gated/3day/';
dataDir_endpoint = '../../data/gated/endpoint/';

numBins = 30;

i = 1;
gene = cell2mat(refList.Gene(i));
baseName_3day = cell2mat(refList.Sample_3day(i));
baseName_endpoint = cell2mat(refList.Sample_Endpoint(i));
[~,label,~] = fileparts(baseName_3day);

%%
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
varNames = {'Red','Green'};
DATA_neg_full = table(Data{1}, Data{2}, 'VariableNames',varNames);
DATA_pos_full = table(Data{3}, Data{4}, 'VariableNames',varNames);

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

% Countour plots
n_grid = 30;
sox2_min = 10;
sox2_max = 10^5;
oct4_min = 10;
oct4_max = 10^5;

figs(1) = figure();
centers = {10.^linspace(log10(oct4_min),log10(oct4_max),n_grid), 10.^linspace(log10(sox2_min),log10(sox2_max),n_grid)};
[n,c] = hist3([DATA_neg.Red, DATA_neg.Green], centers);
contour(c{1},c{2},n',10, 'LineColor', 'black', 'LineWidth', 0.5)
hold on
[n2,c2] = hist3([DATA_pos.Red, DATA_pos.Green], centers);
contour(c2{1},c2{2},n2',10, 'LineColor', [0 0.5 1], 'LineWidth', 1)
legend('CFP-', 'CFP+')
title('Overlaid endpoint of CFP+ and CFP- population')
ylabel('SOX2::YFP (a.u.)','FontSize',10);
xlabel('OCT4::YFP (a.u.)','FontSize',10);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

% Take the log of each data point and add a pseudocount to deal with zeros
pseudocount = 0.00001;
DATA_neg.Red   = log(DATA_neg.Red + pseudocount);
DATA_neg.Green = log(DATA_neg.Green + pseudocount);
DATA_pos.Red   = log(DATA_pos.Red + pseudocount);
DATA_pos.Green = log(DATA_pos.Green + pseudocount);

% Take the ratio of the logs
DATA_neg.Ratio = DATA_neg.Red ./ DATA_neg.Green;
DATA_pos.Ratio = DATA_pos.Red ./ DATA_pos.Green;

% Create histogram data
h = histogram(DATA_neg.Ratio, numBins);
hist_neg = h.Values;
negEdges = h.BinEdges;
h2 = histogram(DATA_pos.Ratio, 'BinEdges', negEdges);
hist_pos = h2.Values;

% Calculate KL divergence
results{3} = f_kldiv(hist_pos, hist_neg, 0.000001);

% Plot histograms
figs(2) = figure();
subplot(1, 2, 1)
histogram(DATA_neg.Ratio, numBins, 'Normalization', 'probability')
title('CFP-negative')
xlabel('OCT4:SOX2')
ylabel('Probability')
subplot(2, 2, 2)
histogram(DATA_pos.Ratio, numBins, 'Normalization', 'probability')
title('CFP-positive')
xlabel('OCT4:SOX2')
ylabel('Probability')

% Endpoint analysis: 

Data = load([dataDir_endpoint baseName_endpoint]);
Data = Data.Data;
    % Cell 1: CFP_negative OCT4
    % Cell 2: CFP_negative SOX2
    % Cell 3: CFP_positive OCT4
    % Cell 4: CFP_positive SOX2
varNames = {'Red','Green'};
DATA_neg_full = table(Data{1}, Data{2}, 'VariableNames',varNames);
DATA_pos_full = table(Data{3}, Data{4}, 'VariableNames',varNames);

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
sox2_max = 10^5;
oct4_min = 10;
oct4_max = 2* 10^4;

figs(3) = figure();
centers = {10.^linspace(log10(oct4_min),log10(oct4_max),n_grid), 10.^linspace(log10(sox2_min),log10(sox2_max),n_grid)};
[n,c] = hist3([DATA_neg.Red, DATA_neg.Green], centers);
contour(c{1},c{2},n',10, 'LineColor', 'black', 'LineWidth', 0.5)
hold on
[n2,c2] = hist3([DATA_pos.Red, DATA_pos.Green], centers);
contour(c2{1},c2{2},n2',10, 'LineColor', [0 0.5 1], 'LineWidth', 1)
title('Overlaid endpoint of CFP+ and CFP- population')
ylabel('SOX2::YFP (a.u.)','FontSize',10);
xlabel('OCT4::YFP (a.u.)','FontSize',10);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

% Gating 

poly_file = ['polys/', label, '_BE.mat'];
if exist(poly_file, 'file')
    load(poly_file, 'coords_BE')
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
    save(poly_file, 'coords_BE')
end
hold off 

% Index BE cells in both CFP+ and CFP- populations
cfp_neg_BE = inpolygon(DATA_neg.Red, DATA_neg.Green, coords_BE(:,1), coords_BE(:,2));
cfp_pos_BE = inpolygon(DATA_pos.Red, DATA_pos.Green, coords_BE(:,1), coords_BE(:,2));

% Calculate BE proportions
BE_neg = sum(cfp_neg_BE) / length(cfp_neg_BE);
BE_pos = sum(cfp_pos_BE) / length(cfp_pos_BE);
% Calculate log of ratios
results{4} = log2( (BE_pos / (1-BE_pos)) / (BE_neg / (1-BE_neg)));

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