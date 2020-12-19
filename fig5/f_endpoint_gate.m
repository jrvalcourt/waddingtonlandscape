function [figs, Data] = f_endpoint_gate(fcsdat, fcshdr, name)
% This function takes in FCS data and returns:
    % The following figures:
        % Figure 1: CFP+/- gating
        % Figure 2: Endpoint gating of CFP- population (scatter + contour)
        % Figure 3: Endpoint gating of CFP+ population (scatter + contour)
    % A data table with raw fluorescence levels of:
        % Column 1: CFP- OCT4::RFP
        % Column 2: CFP- SOX2::YFP
        % Column 3: CFP+ OCT4::RFP
        % Column 4: CFP- SOX2::YFP
    % A cell array with end fate proportions determined by manual gating:
        % Cell 1: Sample name
        % Cell 2: CFP- bipotent ectoderm proportion
        % Cell 3: CFP+ bipotent ectoderm proportion
        % Cell 4: CFP- mesendoderm proportion
        % Cell 5: CFP+ mesendoderm proportion
        
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

% Choose CFP+ and CFP- populations
figs(1) = figure();
scatter(fl_cfp(is_cell), fsc_a(is_cell), '.', 'MarkerFaceAlpha', 0.05,'MarkerEdgeAlpha', 0.05);
hold on
xlabel('mCFP-A','FontSize',20);
ylabel('FSC-A','FontSize',20);
set(gca,'xscale','log')
h = impoly();
coords_cfp_minus = getPosition(h);
h1 = impoly();
coords_cfp_plus = getPosition(h1);
hold off
    
cfp_negative = is_cell & inpolygon(fl_cfp, fsc_a, coords_cfp_minus(:,1), coords_cfp_minus(:,2));
cfp_positive = is_cell & inpolygon(fl_cfp, fsc_a, coords_cfp_plus(:,1), coords_cfp_plus(:,2));

% Output CFP-positive and CFP-negative RFP/YFP values
X1 = fl_oct4(cfp_negative == 1);
Y1 = fl_sox2(cfp_negative == 1);
X2 = fl_oct4(cfp_positive == 1);
Y2 = fl_sox2(cfp_positive == 1);

Data = cell(1,4);
Data{1} = X1;
Data{2} = Y1;
Data{3} = X2;
Data{4} = Y2;

end