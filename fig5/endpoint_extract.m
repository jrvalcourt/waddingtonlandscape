clf; close all; clear variables;
warning('off','all')
warning

dataDir = '../data/mat_files_endpoint/';
dataFiles = dir(fullfile(dataDir, '*.fcs'));
%%

for i=1:length(dataFiles)
    baseName = dataFiles(i).name;
    [~,label,ext] = fileparts(baseName);
    disp([num2str(i), ' out of ', num2str(length(dataFiles)), ', ' label])
    fileName = fullfile(dataDir, baseName);
    [fcsdat, fcshdr, fcsdatscaled, fcsdat_comp] = fca_readfcs(fileName);
    [figs, Data] = f_endpoint_gate(fcsdat, fcshdr, baseName);
    % Data is a cell array of the following order:
        % Cell 1: CFP_negative OCT4
        % Cell 2: CFP_negative SOX2
        % Cell 3: CFP_positive OCT4
        % Cell 4: CFP_positive SOX2
    save(['../data/gated/endpoint/' label '.mat'], 'Data')
    saveas(figs(1), ['../plots/endpoint/gating' label '_cfp_gating.png']);
    close all;
end
