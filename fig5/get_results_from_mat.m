close all
clear all

refList = readtable('reference_list_replicates.xlsx');
dataDir_3day = '../data/mat_files_3day/';
dataDir_endpoint = '../data/mat_files_endpoint/';

numBins = 30;

varNames = {'Gene','Sample', 'KL', 'Log_RR'};
allResults = cell2table(cell(0,4), 'VariableNames', varNames);

fileID = fopen('results_temp.csv','w');

header = 'ii,gene,delta_cutoff,biggest_diff,actual_diff,predicted_diff,sigmoidB_neg,sigmoidB_pos,KL,bias\n';
fprintf(fileID, header);

for i = 1:height(refList)
    gene = cell2mat(refList.Gene(i));
    baseName_3day = cell2mat(refList.Sample_3day(i));
    baseName_endpoint = cell2mat(refList.Sample_Endpoint(i));
    [~, label, ~] = fileparts(baseName_3day);
    [figs, results] = f_results_from_mat(gene, label, dataDir_3day, ...
                                         dataDir_endpoint, baseName_3day, ...
                                         baseName_endpoint, numBins);
    allResults = [allResults;results];
    saveas(figs(1), ['../plots/composite_mat/' label '_3day_contour.png']);
    saveas(figs(2), ['../plots/composite_mat/' label '_3day_histograms.png']);
    saveas(figs(3), ['../plots/composite_mat/' label '_endpoint_contour.png']);
    close all;
    
    
    % do some stuff for jim
    [figs_jim, results_jim] = f_results_from_mat_jim_version2( ...
        gene, label, dataDir_3day, ...
        dataDir_endpoint, baseName_3day, ...
        baseName_endpoint, numBins);
    %allResults = [allResults;results_jim];
    saveas(figs_jim(1), ['../plots/jim_figs/' label '_3day_neg.png']);
    saveas(figs_jim(2), ['../plots/jim_figs/' label '_3day_pos.png']);
    saveas(figs_jim(3), ['../plots/jim_figs/' label '_after_neg.png']);
    saveas(figs_jim(4), ['../plots/jim_figs/' label '_after_pos.png']);
    saveas(figs_jim(5), ['../plots/jim_figs/' label '_3day_neg_lines.png']);
    saveas(figs_jim(6), ['../plots/jim_figs/' label '_3day_pos_lines.png']);
    saveas(figs_jim(7), ['../plots/jim_figs/' label '_ratio_hists.png']);
    saveas(figs_jim(8), ['../plots/jim_figs/' label '_ratio_cdf.png']);
    saveas(figs_jim(9), ['../plots/jim_figs/' label '_ratio_cdf_notlog.png']);
    close all;
    
    diff_results = results_jim{6};
    sigmoid_shifts = results_jim{7};
    
    text_output = strcat(num2str(i), ',', gene, ',', num2str(results_jim{5}), ...
                 ',', num2str(diff_results(1)), ...
                 ',', num2str(diff_results(2)), ...
                 ',', num2str(diff_results(3)), ...
                 ',', num2str(sigmoid_shifts(1)), ...
                 ',', num2str(sigmoid_shifts(2)), ...
                 ',', num2str(results{3}), ...
                 ',', num2str(results{4}), ...
                 '\n');
             
    fprintf(fileID, text_output);
    
    cb_fig = figure();
    colorbar()
    saveas(cb_fig, '../plots/jim_figs/colorbar.png');
    close all;
end

fclose(fileID);

%% Find CFP results

refList = readtable('reference_list_cfp.xlsx');
dataDir_3day = '../data/mat_files_3day/';
dataDir_endpoint = '../data/mat_files_endpoint/';

numBins = 30;

varNames = {'Gene','Sample', 'KL', 'Log_RR'};
cfpResults = cell2table(cell(0,4), 'VariableNames', varNames);

for i = 1:height(refList)
    gene = cell2mat(refList.Gene(i));
    baseName_3day = cell2mat(refList.Sample_3day(i));
    baseName_endpoint = cell2mat(refList.Sample_Endpoint(i));
    [~,label,~] = fileparts(baseName_3day);
    [~, results] = f_results_from_mat(gene, label, dataDir_3day, dataDir_endpoint, baseName_3day, baseName_endpoint, numBins);
    cfpResults = [cfpResults;results];
    close all;
end

kldiv_mean = mean(cfpResults.KL);
kldiv_stdev = std(cfpResults.KL);
logrr_mean = mean(cfpResults.Log_RR);
logrr_stdev = std(cfpResults.Log_RR);


%%

% set up the figure dimensions
fig = figure();
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 8];

% things inside this box should not be labled. probably can set this at 
% the mean CFP +/- 2 sigma
labelx = [logrr_mean - 3*logrr_stdev, logrr_mean + 3*logrr_stdev];
labely = [kldiv_mean - 3*kldiv_stdev, kldiv_mean + 3*kldiv_stdev];

% scatter all points
scatter(allResults.Log_RR, allResults.KL, 'filled')
hold on

% get the labels for each point
sampleNames = allResults.Gene;
dataLabels = cellstr(sampleNames);

% select only the points we want to label
selection = (allResults.Log_RR < labelx(1) | ...
             allResults.Log_RR > labelx(2)) | ...
            (allResults.KL < labely(1) | ...
             allResults.KL > labely(2));
x_to_label = allResults.Log_RR(selection);
y_to_label = allResults.KL(selection);
sampleNames = sampleNames(selection);
dataLabels = dataLabels(selection);

% draw the box
expansion_factor = 1.1;
minx = min(allResults.Log_RR) * expansion_factor;
maxx = max(allResults.Log_RR) * expansion_factor;
miny = min(allResults.KL) * expansion_factor;
maxy = max(allResults.KL) * expansion_factor;
plot([minx, maxx], [labely(2), labely(2)], 'k--', 'LineWidth', 3)
plot([labelx(1), labelx(1)], [miny, maxy], 'k--', 'LineWidth', 3)
plot([labelx(2), labelx(2)], [miny, maxy], 'k--', 'LineWidth', 3)

% finish the plot
set(gca,'linewidth', 2)
xlim([minx, maxx])
ylim([miny, maxy])
xlabel('Log(BE:ME_p_o_s/BE:ME_n_e_g)')
ylabel('KL divergence')

% save without labels
print('../plots/ratio_v_kl_no_labels.png','-dpng','-r600')

% label those pointssssssss
dx = 0.05; dy = 0.05; 
text(x_to_label+dx, y_to_label+dy, dataLabels);

% save with labels
print('../plots/ratio_v_kl.png','-dpng','-r600')
close();