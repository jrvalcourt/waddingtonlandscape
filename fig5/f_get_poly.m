function coords = f_get_poly(label, me_or_ecto, ...
                             data_neg, data_pos, other_poly)

if strcmp(me_or_ecto, 'ecto')
    poly_file = strcat('polys/', label, '_BE.mat');
else
    poly_file = strcat('polys/', label, '_ME.mat');
end

if ~exist(poly_file, 'file')
    
    % params for plot
    n_grid = 30;
    sox2_min = 10;
    sox2_max = 2 .* 10^5;
    oct4_min = 10;
    oct4_max = 10^5;
    
    % get a manageable subset of CFP-negative cells
    if height(data_neg) > 50000
        rand_idx = randsample(height(data_neg), 50000);
        data_neg_disp = data_neg(rand_idx,:);
    else
        data_neg_disp = data_neg;
    end
    
    % plot both cfp+ and cfp- populations
    figure()
    scatter(data_neg_disp.Red, data_neg_disp.Green, '.', ...
        'MarkerEdgeColor','black', 'MarkerFaceAlpha', ...
        0.03,'MarkerEdgeAlpha', 0.03);
    hold on
    scatter(data_pos.Red, data_pos.Green, '.', ...
        'MarkerEdgeColor',[0 0.5 1], 'MarkerFaceAlpha', 0.03, ...
        'MarkerEdgeAlpha', 0.03);
    centers = {10.^linspace(log10(oct4_min),log10(oct4_max),n_grid),...
        10.^linspace(log10(sox2_min),log10(sox2_max),n_grid)};
    [n,c] = hist3([data_neg.Red, data_neg.Green], centers);
    contour(c{1},c{2},n',10, 'LineColor', 'black', 'LineWidth', 0.5)
    [n2,c2] = hist3([data_pos.Red, data_pos.Green], centers);
    contour(c2{1},c2{2},n2',10, 'LineStyle', '--', ...
        'LineColor', [0 0.5 1], 'LineWidth', 1)
    if strcmp(me_or_ecto, 'ecto')
        title('Draw ectoderm polygon:')
    else
        title('Draw mesendoderm polygon:')
    end
    ylabel('SOX2:YFP','FontSize',10);
    xlabel('OCT4:YFP','FontSize',10);
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    
    % draw the (optional) other poly
    if exist('other_poly','var')
        plot(other_poly(:,1), other_poly(:,2))
    end
    
    % let the user draw the polygon
    h = drawpolygon();
    coords = h.Position;
    save(poly_file, 'coords')
    
    close()
end
load(poly_file, 'coords')
end