clf; close all; clear variables;
warning('off','all')
warning

dataDir = '../data/mat_files_es/';
dataFiles = dir(fullfile(dataDir, '*.fcs'));

%% All data is plotted in logspace but kept in linear format 

%for i = 2:length(dataFiles)
for i = 4:5
    baseName = dataFiles(i).name;
    [~,label,ext] = fileparts(baseName);
    fileName = fullfile(dataDir, baseName);
    [fcsdat, fcshdr, fcsdatscaled, fcsdat_comp] = fca_readfcs(fileName);
    [figs, Data] = f_es_gate(fcsdat, fcshdr, baseName);
    save(['../data/gated/3day_es/' label '.mat'], 'Data')
    saveas(figs(1), ['../plots/3day_es/' label '_es.png']);
    close all;
end
