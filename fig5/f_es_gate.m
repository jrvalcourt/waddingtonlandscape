function [figs, Data] = f_es_gate(fcsdat, fcshdr, name)
% This function takes in FCS data and returns:
    % The following figures:
        % Figure 1: CFP+/- gating
        % Figure 2: Contour subplots of CFP- (top) and CFP+ (bottom) 
    % A data table with raw fluorescence levels of:
        % Column 1: CFP- OCT4::RFP
        % Column 2: CFP- SOX2::YFP
        % Column 3: CFP+ OCT4::RFP
        % Column 4: CFP- SOX2::YFP

figs = [];

% Pull out the relevant channels.
index_fsc_a = find(strcmp({fcshdr.par.name}, 'FSC-A')==1);
index_fsc_h = find(strcmp({fcshdr.par.name}, 'FSC-H')==1);
index_fsc_w = find(strcmp({fcshdr.par.name}, 'FSC-W')==1);
index_ssc_a = find(strcmp({fcshdr.par.name}, 'SSC-A')==1);
index_ssc_h = find(strcmp({fcshdr.par.name}, 'SSC-H')==1);
index_ssc_w = find(strcmp({fcshdr.par.name}, 'SSC-W')==1);
index_sox2 = find(strcmp({fcshdr.par.name}, 'FITC-A')==1);
index_oct4 = find(strcmp({fcshdr.par.name}, 'PE-Texas Red-A')==1);
index_cfp = find(strcmp({fcshdr.par.name}, 'mCFP-A')==1);
fsc_a = fcsdat(:,index_fsc_a);
fsc_h = fcsdat(:,index_fsc_h);
fsc_w = fcsdat(:,index_fsc_w);
ssc_a = fcsdat(:,index_ssc_a);
ssc_h = fcsdat(:,index_ssc_h);
ssc_w = fcsdat(:,index_ssc_w);
fl_sox2 = fcsdat(:,index_sox2);
fl_oct4 = fcsdat(:,index_oct4);
fl_cfp = fcsdat(:,index_cfp);

is_valid = (fsc_a < 2 * 10.^5) & (ssc_a < 2 * 10.^5);

% Choose cells
fsc_ssc_poly_file = ['polys/', name, '_fsc_ssc.mat'];
if exist(fsc_ssc_poly_file, 'file')
    load(fsc_ssc_poly_file, 'coords1')
else
    %view = is_valid(1:10000);
    scatter(fsc_a(is_valid), ssc_a(is_valid), '.', 'MarkerFaceAlpha', 0.05,'MarkerEdgeAlpha', 0.05);
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

% Choose single cells
fscw_fsch_poly_file = ['polys/', name, '_fscw_fsch.mat'];
if exist(fscw_fsch_poly_file, 'file')
    load(fscw_fsch_poly_file, 'coords2')
else
    hold off
    scatter(fsc_w(is_right_fsc_ssc), fsc_h(is_right_fsc_ssc), '.', 'MarkerFaceAlpha', 0.05,'MarkerEdgeAlpha', 0.05);
    hold on
    xlabel('FSC-W','FontSize',20);
    ylabel('FSC-H','FontSize',20);
    h = impoly();
    coords2 = getPosition(h);
    save(fscw_fsch_poly_file, 'coords2')
end
hold off 

is_right_fscw_fsch = inpolygon(fsc_w, fsc_h, coords2(:,1), coords2(:,2));
is_cell = is_right_fsc_ssc & is_right_fscw_fsch;

% Gate cells that are not SOX2-silenced
figure()
scatter(fl_oct4(is_cell), fl_sox2(is_cell), '.', 'MarkerFaceAlpha', 0.05,'MarkerEdgeAlpha', 0.05);
hold on
set(gca,'xscale','log')
set(gca,'yscale','log')
title('Single cells')
ylabel('SOX2::YFP','FontSize',10);
xlabel('OCT4::RFP','FontSize',10);
h = impoly();
coords_cells_want = getPosition(h);
hold off

is_cell = is_cell & inpolygon(fl_oct4, fl_sox2, coords_cells_want(:,1), coords_cells_want(:,2));

% For plotting contours
n_grid = 30;
sox2_min = 10;
sox2_max = 1 * 10^5;
oct4_min = 4 * 10^2;
oct4_max = 2 * 10^5;

% Final contour
figs(1) = figure();
%scatter(fl_oct4(cfp_negative), fl_sox2(cfp_negative), '.', 'MarkerFaceAlpha', 0.05,'MarkerEdgeAlpha', 0.05);
centers = {10.^linspace(log10(sox2_min),log10(sox2_max),n_grid), 10.^linspace(log10(oct4_min),log10(oct4_max),n_grid)};
[n,c] = hist3([fl_oct4(is_cell), fl_sox2(is_cell)], centers);
contour(c{1},c{2},n',10, 'LineWidth', 1)
hold on
set(gca,'xscale','log')
set(gca,'yscale','log')
title('ES cells for calibration')
ylabel('SOX2::YFP','FontSize',10);
xlabel('OCT4::RFP','FontSize',10);
hold off

% Output CFP-positive and CFP-negative RFP/YFP tables
X1 = fl_oct4(is_cell == 1);
Y1 = fl_sox2(is_cell == 1);

Data = cell(1,2);
Data{1} = X1;
Data{2} = Y1;

end