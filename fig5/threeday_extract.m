clf; close all; clear variables;
warning('off','all')
warning

dataDir = '../data/mat_files_3day/';
dataFiles = dir(fullfile(dataDir, '*.fcs'));

%% All data is plotted in logspace but kept in linear format 

%for i = 1:length(dataFiles)
for i = 1
    baseName = dataFiles(i).name;
    [~,label,ext] = fileparts(baseName);
    fileName = fullfile(dataDir, baseName);
    [fcsdat, fcshdr, fcsdatscaled, fcsdat_comp] = fca_readfcs(fileName);
    [figs, Data] = f_3day_gate(fcsdat, fcshdr, baseName);
    save(['../data/gated/3day/' label '.mat'], 'Data')
    saveas(figs(1), ['../plots/' label '_cfp.png']);
    saveas(figs(2), ['../plots/' label '_oct4_sox2.png']);
    close all;
end
