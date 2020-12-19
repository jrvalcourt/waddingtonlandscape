%https://www.mathworks.com/matlabcentral/answers/194554-how-can-i-use-and-display-two-different-colormaps-on-the-same-figure
%https://www.mathworks.com/help/matlab/ref/colormapeditor.html

clf; close all; clear all;

rng(2);

sigmoidA = 8.37755655;

% load custom colormaps
load('colormaps/purples.mat');
load('colormaps/greens.mat');
load('colormaps/yellows.mat');
load('colormaps/blues.mat');
cmaps = {greens, purples, blues, yellows};
purple_color = [102./255, 45./255, 145./255];
green_color  = [16./255, 167./255,  74./255];
blue_color   = [40./255, 126./255, 194./255];

% create containers
c = containers.Map;
channel_set = containers.Map;

[p4_x, p4_y, p7_x, p7_y] = get_gates(['data/EXP0146.xml']);

% file paths
c('EXP0252_ES')    = 'data/EXP0252_mSR2_pre-bmp-acta_es_accutase.fcs';
c('EXP0252_3days') = 'data/EXP0252_mSR2_pre-bmp-acta_no-transduction.fcs';
c('EXP0252_after') = 'data/EXP0252_mSR2_after-bmp-acta_no-transduction.fcs';
c('EXP0249_ES')    = 'data/EXP0249_mSR2_pre-bmp-acta_es_accutase.fcs';
c('EXP0249_3days') = 'data/EXP0249_mSR2_pre-bmp-acta_no_virus.fcs';
c('EXP0249_after') = 'data/EXP0249_mSR2_after-bmp-acta_no_virus.fcs';
c('EXP0146_3days') = 'data/EXP0146_mSR2_ChIP-seq_sample.fcs';

% channels to use
% fsc_a fsc_h fsc_w ssc_a fl_sox2 fl_oct4
channel_set('EXP0249_ES')    = [1,2,3,4,9,14];
channel_set('EXP0146_3days') = [2,3,4,5,8,9];
channel_set('EXP0249_after') = [1,2,3,4,9,14];
channel_set('EXP0249_3days') = [1,2,3,4,9,14];
channel_set('EXP0252_ES')    = [1,2,3,4,9,14];
channel_set('EXP0252_3days') = [1,2,3,4,9,14];
channel_set('EXP0252_after') = [1,2,3,4,9,14];

es_key        = 'EXP0249_ES';
d3_key        = 'EXP0249_3days';
after_key     = 'EXP0249_after';
sorted_key    = 'EXP0146_3days';
es_alt_key    = 'EXP0252_ES';
d3_alt_key    = 'EXP0252_3days';
after_alt_key = 'EXP0252_after';

% get data for the plot with 3 populations
[es_oct4,    es_sox2]      = get_data(get_fcsdat(c(es_key)), ...
    es_key, channel_set(es_key));
[sorted_oct4, sorted_sox2] = get_data(get_fcsdat(c(sorted_key)), ...
    sorted_key, channel_set(sorted_key));
[d3_oct4,    d3_sox2]      = get_data(get_fcsdat(c(d3_key)), ...
    d3_key, channel_set(d3_key));
[after_oct4, after_sox2]   = get_data(get_fcsdat(c(after_key)), ...
    after_key, channel_set(after_key));

% alt populations
[es_oct4_alt, es_sox2_alt] = get_data(get_fcsdat(c(es_alt_key)), ...
    es_alt_key, channel_set(es_alt_key));
[d3_oct4_alt, d3_sox2_alt] = get_data(get_fcsdat(c(d3_alt_key)), ...
    d3_alt_key, channel_set(d3_alt_key));
[after_oct4_alt, after_sox2_alt] = get_data(get_fcsdat(c(after_alt_key)), ...
    after_alt_key, channel_set(after_alt_key));

% trim silenced cells -- in some cells, the SOX2 reporter gets silenced,
% which presents as a small population with extremely low SOX2:YFP signal.
% we do not show these silenced cells in this figure for visual clarity.
not_silenced = d3_sox2_alt > 2^13.5;
d3_sox2_alt = d3_sox2_alt(not_silenced);
d3_oct4_alt = d3_oct4_alt(not_silenced);

not_silenced = sorted_sox2 > 2^12.5;
sorted_sox2 = sorted_sox2(not_silenced);
sorted_oct4 = sorted_oct4(not_silenced);

% do gaussian mixture
[ecto, mesendo, frac_me] = gaussianmixture(after_oct4, after_sox2);
[ecto_alt, mesendo_alt, other1_alt, other2_alt, frac_me_alt] = ...
    gaussianmixture4(after_oct4_alt, after_sox2_alt);

% plot the sorted population
mplot = plot_colored_contours_no_fill({sorted_oct4}, ...
    {sorted_sox2}, {purples});
hold on;
plot(10.^p4_x, 10.^p4_y, 'k', 'LineWidth', 3)
plot(10.^p7_x, 10.^p7_y, 'k', 'LineWidth', 3)
saveas(mplot, 'plots/sorted_pop_day3_purple.png')

% make plots with the 4 populations
mplot = plot_four_contours({es_oct4, d3_oct4, ecto(:,1), mesendo(:,1)}, ...
    {es_sox2, d3_sox2, ecto(:,2), mesendo(:,2)}, ...
    cmaps);
saveas(mplot, 'plots/es-and-d3_oct4-sox2_sep-me-biecto.png')

% norm the values to es
norm_es_oct4 = es_oct4 ./ mean(es_oct4);
norm_es_sox2 = es_sox2 ./ mean(es_sox2);
norm_d3_oct4 = d3_oct4 ./ mean(es_oct4);
norm_d3_sox2 = d3_sox2 ./ mean(es_sox2);
norm_ecto_oct4 = ecto(:,1) ./ mean(es_oct4);
norm_ecto_sox2 = ecto(:,2) ./ mean(es_sox2);

% norm the values to es
norm_es_oct4_alt = es_oct4_alt ./ mean(es_oct4_alt);
norm_es_sox2_alt = es_sox2_alt ./ mean(es_sox2_alt);
norm_d3_oct4_alt = d3_oct4_alt ./ mean(es_oct4_alt);
norm_d3_sox2_alt = d3_sox2_alt ./ smean(es_sox2_alt);
norm_ecto_oct4_alt = ecto_alt(:,1) ./ mean(es_oct4_alt);
norm_ecto_sox2_alt = ecto_alt(:,2) ./ mean(es_sox2_alt);

% find the sigmoid location
log2sigmoidloc = f_infer_sigmoid_loc(log2(norm_d3_oct4 ./ norm_d3_sox2), ...
    frac_me, sigmoidA);
log2sigmoidloc_alt = f_infer_sigmoid_loc(...
    log2(norm_d3_oct4_alt ./ norm_d3_sox2_alt), frac_me_alt, sigmoidA);

% make some plots
mplot = plot_three_hists({norm_es_oct4_alt, norm_d3_oct4_alt, norm_ecto_oct4_alt}, ...
                         {norm_es_sox2_alt, norm_d3_sox2_alt, norm_ecto_sox2_alt}, ...
                         {green_color, purple_color, blue_color}, ...
                         sigmoidA, log2sigmoidloc_alt);
saveas(mplot, 'plots/three_hists.pdf')    

mplot = plot_hist_with_sigmoid(norm_d3_oct4_alt, norm_d3_sox2_alt, ...
                         purple_color, sigmoidA, log2sigmoidloc_alt);
saveas(mplot, 'plots/hist_with_sigmoid.pdf')    

mplot = plot_hist_with_gates(sorted_oct4, sorted_sox2, purple_color, ...
                         10.^p4_x, 10.^p4_y, 10.^p7_x, 10.^p7_y);
saveas(mplot, 'plots/hist_with_gates.pdf')  

close all;

% functionssssss
function [p4_x, p4_y, p7_x, p7_y] = get_gates(fname)
expt_data = parseXML(fname);

p4_points = expt_data.Children.Children(9).Children.Children(16).Children(5).Children(8).Children.Children;
p7_points = expt_data.Children.Children(9).Children.Children(16).Children(8).Children(8).Children.Children;

p4_x = [];
p4_y = [];
for ii=1:length(p4_points)
    p4_x = [p4_x, str2double(p4_points(ii).Attributes(1).Value)];
    p4_y = [p4_y, str2double(p4_points(ii).Attributes(2).Value)];
end
p4_x = [p4_x, str2double(p4_points(1).Attributes(1).Value)];
p4_y = [p4_y, str2double(p4_points(1).Attributes(2).Value)];

p7_x = [];
p7_y = [];
for ii=1:length(p7_points)
    p7_x = [p7_x, str2double(p7_points(ii).Attributes(1).Value)];
    p7_y = [p7_y, str2double(p7_points(ii).Attributes(2).Value)];
end
p7_x = [p7_x, str2double(p7_points(1).Attributes(1).Value)];
p7_y = [p7_y, str2double(p7_points(1).Attributes(2).Value)];
end
function fcsdat = get_fcsdat(fname)
[fcsdat, ~, ~, ~] = fca_readfcs(fname);
%[fcsdat, fcshdr, fcsdatscaled, fcsdat_comp] = fca_readfcs(fname);
%x = 1+1;
end
function fig = plot_colored_contours_no_fill(oct4s, sox2s, cmaps)
n_grid = 30;
n_levels = 10;
amap = 'rampup';
lw = 2;
xlims = [1.5 * 10^2, 3 * 10^4];
ylims = [3 * 10^3, 9 * 10^4];

oct4_min = 1 * 10^0;
oct4_max = 3 * 10^5;

sox2_min = 1 * 10^0;
sox2_max = 3 * 10^5;

fig = figure();
centers = {10.^linspace(log10(sox2_min),log10(sox2_max),n_grid), ...
    10.^linspace(log10(oct4_min),log10(oct4_max),n_grid)};

hold on;
contours = {};

ax1 = axes;
fl_oct4 = oct4s{1};
fl_sox2 = sox2s{1};
[n,c] = hist3([fl_oct4, fl_sox2], centers);
% https://www.mathworks.com/matlabcentral/answers/250279-how-to-make-contour-plots-transparent-in-matlab-r2015a
[~, contours{1}] = contour(c{1}, c{2}, n', n_levels, 'LineWidth', lw);
%set(contours{1},'LineColor','none')
xlim(xlims)
ylim(ylims)
set(ax1,'xscale','log')
set(ax1,'yscale','log')
colormap(ax1, cmaps{1})

ax1.Visible = 'off';

set(ax1,'Position',[.17 .11 .685 .79]);
cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
end
function fig = plot_four_contours(oct4s, sox2s, cmaps)
n_grid = 60;
n_levels = 20;
amap = 'rampup';
lw = 1.5;
xlims = [1 * 10^1, 7 * 10^4];
ylims = [3 * 10^2, 4 * 10^4];

oct4_min = 1 * 10^0;
oct4_max = 3 * 10^5;

sox2_min = 1 * 10^0;
sox2_max = 3 * 10^5;

fig = figure();
centers = {10.^linspace(log10(sox2_min),log10(sox2_max),n_grid), ...
    10.^linspace(log10(oct4_min),log10(oct4_max),n_grid)};

hold on;
contours = {};

ax1 = axes;
fl_oct4 = oct4s{1};
fl_sox2 = sox2s{1};
[n,c] = hist3([fl_oct4, fl_sox2], centers);
% https://www.mathworks.com/matlabcentral/answers/250279-how-to-make-contour-plots-transparent-in-matlab-r2015a
[~, contours{1}] = contourf(c{1}, c{2}, n', n_levels, 'LineWidth', lw);
set(contours{1},'LineColor','none')
xlim(xlims)
ylim(ylims)
set(ax1,'xscale','log')
set(ax1,'yscale','log')
colormap(ax1, cmaps{1})

ax2 = axes;
fl_oct4 = oct4s{2};
fl_sox2 = sox2s{2};
[n,c] = hist3([fl_oct4, fl_sox2], centers);
% https://www.mathworks.com/matlabcentral/answers/250279-how-to-make-contour-plots-transparent-in-matlab-r2015a
[~, contours{2}] = contourf(c{1}, c{2}, n', n_levels, 'LineWidth', lw);
set(contours{2},'LineColor','none')
xlim(xlims)
ylim(ylims)
set(ax2,'xscale','log')
set(ax2,'yscale','log')
colormap(ax2, cmaps{2})

ax3 = axes;
fl_oct4 = oct4s{3};
fl_sox2 = sox2s{3};
[n,c] = hist3([fl_oct4, fl_sox2], centers);
% https://www.mathworks.com/matlabcentral/answers/250279-how-to-make-contour-plots-transparent-in-matlab-r2015a
[~, contours{3}] = contourf(c{1}, c{2}, n', n_levels, 'LineWidth', lw);
set(contours{3},'LineColor','none')
xlim(xlims)
ylim(ylims)
set(ax3,'xscale','log')
set(ax3,'yscale','log')
colormap(ax3, cmaps{3})

ax4 = axes;
fl_oct4 = oct4s{4};
fl_sox2 = sox2s{4};
[n,c] = hist3([fl_oct4, fl_sox2], centers);
% https://www.mathworks.com/matlabcentral/answers/250279-how-to-make-contour-plots-transparent-in-matlab-r2015a
[~, contours{4}] = contourf(c{1}, c{2}, n', n_levels, 'LineWidth', lw);
set(contours{4},'LineColor','none')
xlim(xlims)
ylim(ylims)
set(ax4,'xscale','log')
set(ax4,'yscale','log')
colormap(ax4, cmaps{4})

linkaxes([ax1, ax2, ax3, ax4]);

%ax1.Visible = 'off';
%ax1.XTick = [];
%ax1.YTick = [];
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];
ax4.Visible = 'off';
ax4.XTick = [];
ax4.YTick = [];

ax4.CLim = ax3.CLim;

set([ax1, ax2, ax3, ax4],'Position',[.17 .11 .685 .79]);
cb1 = colorbar(ax1,'Position',[0.80 0.11 0.05 0.85], 'Ticks',[], 'TickLabels',{});
cb2 = colorbar(ax2,'Position',[0.95 0.11 0.05 0.85], 'Ticks',[], 'TickLabels',{});
cb3 = colorbar(ax3,'Position',[0.90 0.11 0.05 0.85], 'Ticks',[], 'TickLabels',{});
cb4 = colorbar(ax4,'Position',[0.85 0.11 0.05 0.85], 'Ticks',[], 'TickLabels',{});

drawnow;
for ii = 1:numel(contours)
    hContour = contours{ii};
    hFills = hContour.FacePrims;
    [hFills.ColorType] = deal('truecoloralpha');
    for idx = 1 : numel(hFills)
        hFills(idx).ColorData(4) = (idx-1) * (255 / n_levels);   % default=255
    end
end
end
function [oct4, sox2] = get_data(fcsdat, name, channels)

% pull out the relevant channels
fsc_a = fcsdat(:,channels(1));
fsc_h = fcsdat(:,channels(2));
fsc_w = fcsdat(:,channels(3));
ssc_a = fcsdat(:,channels(4));
fl_sox2 = fcsdat(:,channels(5));
fl_oct4 = fcsdat(:,channels(6));

fsc_ssc_poly_file = ['polys/', name, '_fsc_ssc.mat'];
if exist(fsc_ssc_poly_file, 'file')
    load(fsc_ssc_poly_file, 'coords1')
else
    scatter(fsc_a, ssc_a, '.', ...
        'MarkerFaceAlpha', 0.05,'MarkerEdgeAlpha', 0.05);
    hold on
    xlabel('FSC-A','FontSize',20);
    ylabel('SSC-A','FontSize',20);
    xlim([0, 2 * 10.^5])
    ylim([0, 2 * 10.^5])
    h = impoly();
    coords1 = getPosition(h);
    save(fsc_ssc_poly_file, 'coords1')
end

is_right_fsc_ssc = inpolygon(fsc_a, ssc_a, coords1(:,1), coords1(:,2));

fscw_fsch_poly_file = ['polys/', name, '_fscw_fsch.mat'];
if exist(fscw_fsch_poly_file, 'file')
    load(fscw_fsch_poly_file, 'coords2')
else
    hold off
    scatter(fsc_w(is_right_fsc_ssc), fsc_h(is_right_fsc_ssc), ...
        '.', 'MarkerFaceAlpha', 0.05,'MarkerEdgeAlpha', 0.05);
    hold on
    xlabel('FSC-W','FontSize',20);
    ylabel('FSC-H','FontSize',20);
    h = impoly();
    coords2 = getPosition(h);
    save(fscw_fsch_poly_file, 'coords2')
end

is_right_fscw_fsch = inpolygon(fsc_w, fsc_h, coords2(:,1), coords2(:,2));

is_cell = is_right_fsc_ssc & is_right_fscw_fsch;

oct4 = fl_oct4(is_cell);
sox2 = fl_sox2(is_cell);

nonzero = (oct4 > 0) & (sox2 > 0);
oct4 = oct4(nonzero);
sox2 = sox2(nonzero);
end
function fig = plot_three_hists(oct4s, sox2s, colors, sigmoidA, sigmoidB)

nbins=100;

fig = figure('Renderer', 'painters', 'Position', [10 10 600 300]);
es_log2ratio   = log2(oct4s{1} ./ sox2s{1});
d3_log2ratio   = log2(oct4s{2} ./ sox2s{2});
ecto_log2ratio = log2(oct4s{3} ./ sox2s{3});

xs = linspace(min(ecto_log2ratio), max(es_log2ratio), nbins);
xcenters = xs(1:(length(xs)-1)) + 0.5 .* diff(xs);

count_es   = histcounts(es_log2ratio,   xs);
count_d3   = histcounts(d3_log2ratio,   xs);
count_ecto = histcounts(ecto_log2ratio, xs);

prob_es   = count_es   ./ sum(count_es);
prob_d3   = count_d3   ./ sum(count_d3);
prob_ecto = count_ecto ./ sum(count_ecto);

plot(xcenters, prob_es,   'o-', 'Color', colors{1}, ...
    'LineWidth', 3)
hold on
plot(xcenters, prob_d3,   'o-', 'Color', colors{2}, ...
    'LineWidth', 3)
plot(xcenters, prob_ecto, 'o-', 'Color', colors{3}, ...
    'LineWidth', 3)
set(gca,'Visible','off')
xlim([-9,2])
end
function [ecto, mesendo, frac_me] = ...
    gaussianmixture(oct4s, sox2s)
% gaussian mixture
gm = fitgmdist([log2(oct4s), log2(sox2s)], 2);
clusters = cluster(gm, [log2(oct4s), log2(sox2s)]);
cl1 = [oct4s(clusters == 1), sox2s(clusters == 1)];
cl2 = [oct4s(clusters == 2), sox2s(clusters == 2)];
cl1_mean = mean(cl1);
cl2_mean = mean(cl2);
if cl1_mean(1) > cl2_mean(1)
    cl2_temp = cl2;
    cl2 = cl1;
    cl1 = cl2_temp;
end
ecto = cl1;
mesendo = cl2;
frac_me = length(mesendo) ./ (length(ecto) + length(mesendo));
end
function [ecto, mesendo, other1, other2, frac_me] = ...
    gaussianmixture4(oct4s, sox2s)
% gaussian mixture
gm = fitgmdist([log2(oct4s), log2(sox2s)], 4);
clusters = cluster(gm, [log2(oct4s), log2(sox2s)]);
cl1 = [oct4s(clusters == 1), sox2s(clusters == 1)];
cl2 = [oct4s(clusters == 2), sox2s(clusters == 2)];
cl3 = [oct4s(clusters == 3), sox2s(clusters == 3)];
cl4 = [oct4s(clusters == 4), sox2s(clusters == 4)];
cl1_mean = mean(log2(cl1(:,1) ./ cl1(:,2)));
cl2_mean = mean(log2(cl2(:,1) ./ cl2(:,2)));
cl3_mean = mean(log2(cl3(:,1) ./ cl3(:,2)));
cl4_mean = mean(log2(cl4(:,1) ./ cl4(:,2)));

[~, idx] = sort([cl1_mean, cl2_mean, cl3_mean, cl4_mean]);

ecto = [oct4s(clusters == idx(1)), sox2s(clusters == idx(1))];
mesendo = [oct4s(clusters == idx(4)), sox2s(clusters == idx(4))];
other1 = [oct4s(clusters == idx(2)), sox2s(clusters == idx(2))];
other2 = [oct4s(clusters == idx(3)), sox2s(clusters == idx(3))];
frac_me = length(mesendo) ./ (length(ecto) + length(mesendo));
end
function fig = plot_hist_with_sigmoid(oct4, sox2, color, ...
    sigmoidA, sigmoidB)

nbins=40;
xlims = [-8,3];

fig = figure('Renderer', 'painters', 'Position', [10 10 600 300]);
log2ratio = log2(oct4 ./ sox2);
lower = quantile(log2ratio, 0.15);
upper = quantile(log2ratio, 0.85);

xs = linspace(min(log2ratio), max(log2ratio), nbins);
xcenters = xs(1:(length(xs)-1)) + 0.5 .* diff(xs);

counts = histcounts(log2ratio, xs);
probs = counts ./ sum(counts);

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
plot(xcenters, probs, 'o-', 'Color', color, ...
    'LineWidth', 2)
hold on;
xlim(xlims)
set(gca,'Visible','off')

ax2 = axes('Position', ax1_pos, 'Color', 'none', 'YAxisLocation', 'right');
%plot([lower, lower], [0, 1], 'k', 'LineWidth', 2)
hold on;
%plot([upper, upper], [0, 1], 'k', 'LineWidth', 2)
plot(xcenters, f_sigmoid(xcenters, sigmoidA, sigmoidB), ...
    'ko-', 'LineWidth', 2)
xlim(xlims)
%set(gca,'Visible','off')
end
function fig = plot_hist_with_gates(oct4, sox2, color, ...
    p4_x, p4_y, p7_x, p7_y)

nbins=40;
xlims = [-8,0];

fig = figure();
log2ratio = log2(oct4 ./ sox2);

xs = linspace(min(log2ratio), max(log2ratio), nbins);
xcenters = xs(1:(length(xs)-1)) + 0.5 .* diff(xs);

counts = histcounts(log2ratio, xs);
probs = counts ./ sum(counts);

plot(xcenters, probs, 'o-', 'Color', color, ...
    'LineWidth', 2)
hold on;
lower = log2(p4_x(1) ./ p4_y(1));
upper = log2(p7_x(1) ./ p7_y(1));
plot([lower, lower], [0, 0.2], 'k', 'LineWidth', 2)
plot([upper, upper], [0, 0.2], 'k', 'LineWidth', 2)
xlim(xlims)
%set(gca,'Visible','off')

end