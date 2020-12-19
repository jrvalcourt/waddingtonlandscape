close all;
clear all;

rng(1)

sigmoidAs = [1,2,3,4,5,6,7,8, 8.37755655, 9]

for aa=1:length(sigmoidAs)

% params
sigmoidA = sigmoidAs(aa);
nbins = 20;
pseudocount = 0.00001;

% input file locations
refList = readtable('reference_list_full.xlsx');
dataDir_3day  = '../data/mat_files_3day/';
dataDir_endpt = '../data/mat_files_endpoint/';
dataDir_es    = '../data/mat_files_es/';

% output file locations
fileID = fopen('results.csv','w');

% print for 
header = strcat('ii,gene,sample_name,sigmoidB_neg,sigmoidB_pos,', ...
                'frac_me_neg,frac_me_pos,KL,KL_sign,', ...
                'exp_frac_me_without_sigmoid_shift\n');
fprintf(fileID, header);

es_registered_3day_neg_logratios = {};
es_registered_sigmoid_locations = [];
es_registered_es_logratios = {};
sample_names = {};
gene_shifts = containers.Map;
gene2idx = containers.Map;
osr_distributions_neg = {};
osr_distributions_pos = {};

count = 0;
for ii = 1:height(refList)
    
    if refList.skip(ii) > 0.5
        continue
    end
    count = count + 1; 
    
    % parse data for each line
    gene = cell2mat(refList.gene(ii));
    if ~isKey(gene2idx, gene)
        gene2idx(gene) = [];
    end
    gene2idx(gene) = [gene2idx(gene) ii];
    baseName_3day  = cell2mat(refList.sample_3day(ii));
    baseName_endpt = cell2mat(refList.sample_endpoint(ii));
    baseName_es    = cell2mat(refList.matched_es(ii));
    [~, label, ~] = fileparts(baseName_3day);
    parts = split(label, '_');
    temp = strcat(gene, '_', parts(1), '_', parts(2));
    sample_name = temp{1};
    sample_names{count} = sample_name;
    
    % build file paths
    fname_3day  = strcat(dataDir_3day,  baseName_3day);
    fname_endpt = strcat(dataDir_endpt, baseName_endpt);
    fname_es    = strcat(dataDir_es,    baseName_es);

    % get the raw results for this sample
    [data_3day_neg,  data_3day_pos]  = f_fetch_data(fname_3day);
    [data_endpt_neg, data_endpt_pos] = f_fetch_data(fname_endpt);
    
    % gating
    coords_me   = f_get_poly(sample_name, 'mesendo', ...
                             data_endpt_neg, data_endpt_pos);
    coords_ecto = f_get_poly(sample_name, 'ecto', ...
                             data_endpt_neg, data_endpt_pos, coords_me);
    
    % calc fractions
    frac_me_neg = f_get_frac_me(data_endpt_neg, coords_me, coords_ecto);
    frac_me_pos = f_get_frac_me(data_endpt_pos, coords_me, coords_ecto);
    
    if isfile(fname_es)
        
        % get the ES data
        data_es = f_fetch_data_es(fname_es);
        
        % get the ES means for OCT4 and SOX2
        mean_oct4_es = mean(data_es.Red);
        mean_sox2_es = mean(data_es.Green);
        
        % normalize by ES data
        es_registered_es_logratio = ...
            log2(data_es.Red   ./ mean_oct4_es)  - ...
            log2(data_es.Green ./ mean_sox2_es);  
        es_registered_3day_neg_logratio = ...
            log2(data_3day_neg.Red   ./ mean_oct4_es)  - ...
            log2(data_3day_neg.Green ./ mean_sox2_es);
        es_registered_3day_pos_logratio = ...
            log2(data_3day_pos.Red   ./ mean_oct4_es)  - ...
            log2(data_3day_pos.Green ./ mean_sox2_es);
        
        % normalize by ES data
        es_registered_es_ratio = ...
            (data_es.Red   ./ mean_oct4_es) ./ ...
            (data_es.Green ./ mean_sox2_es);  
        es_registered_3day_neg_ratio = ...
            (data_3day_neg.Red   ./ mean_oct4_es) ./ ...
            (data_3day_neg.Green ./ mean_sox2_es);
        es_registered_3day_pos_ratio = ...
            (data_3day_pos.Red   ./ mean_oct4_es) ./ ...
            (data_3day_pos.Green ./ mean_sox2_es);
        
        % infer where the sigmoid is in the normalized coords
        sigmoidB_neg_es_registered = ...
            f_infer_sigmoid_loc(es_registered_3day_neg_logratio, ...
            frac_me_neg, sigmoidA);
        sigmoidB_pos_es_registered = ...
            f_infer_sigmoid_loc(es_registered_3day_pos_logratio, ...
            frac_me_pos, sigmoidA);
        exp_frac_me_from_pOS = f_calc_exp_frac_me( ...
            es_registered_3day_pos_logratio, sigmoidA, ...
            sigmoidB_neg_es_registered);
        
        % plotsssssss
        f_plot_sigmoid_shift(strcat(sample_name, '_esnorm'), sigmoidA, ...
            sigmoidB_neg_es_registered, [sigmoidB_pos_es_registered]);
        f_plot_logratio_hists(strcat(sample_name, '_esnorm'), ...
            log2(es_registered_3day_neg_ratio), ...
            log2(es_registered_3day_pos_ratio));
    end
    
    % calculate log ratios
    logratios_3day_neg = log2(data_3day_neg.Red ./ data_3day_neg.Green);
    logratios_3day_pos = log2(data_3day_pos.Red ./ data_3day_pos.Green);
   
    % infer where the sigmoid is
    sigmoidB_neg = ...
        f_infer_sigmoid_loc(logratios_3day_neg, frac_me_neg, sigmoidA);
    sigmoidB_pos = ...
        f_infer_sigmoid_loc(logratios_3day_pos, frac_me_pos, sigmoidA);
    if ~isKey(gene_shifts, gene)
        gene_shifts(gene) = [];
    end
    gene_shifts(gene) = [gene_shifts(gene), sigmoidB_pos - sigmoidB_neg];
    
    % store the OSR distributions
    osr_distributions_neg{ii} = logratios_3day_neg;
    osr_distributions_pos{ii} = logratios_3day_pos;
    
    % calculate the KL divergence between the cfp+ and cfp- distributions
    % at 3 days
    kl_3day = f_kl_div(logratios_3day_neg, logratios_3day_pos, ...
                       nbins, pseudocount);
    kl_sign = NaN;
    if mean(logratios_3day_pos) > mean(logratios_3day_neg)
        kl_sign = 1;
    else
        kl_sign = -1;
    end
                   
    % plots
%     f_plot_sigmoid_shift(sample_name, sigmoidA, ...
%             sigmoidB_neg, [sigmoidB_pos]);
%     f_plot_logratio_hists(sample_name, ...
%             logratios_3day_neg, ...
%             logratios_3day_pos);
    
    % logratio hists
%     mean_logratio_es = NaN;
%     if ~isfile(fname_es)
%         f_plot_logratio_hists(sample_name, ...
%             logratios_3day_neg, ...
%             logratios_3day_pos);
%         f_plot_prob_colormap(strcat(sample_name, '_neg'), ...
%             sigmoidA, sigmoidB_neg);
%         f_plot_prob_colormap(strcat(sample_name, '_pos'), ...
%             sigmoidA, sigmoidB_pos);
%         f_plot_sigmoid_shift(sample_name, sigmoidA, ...
%             sigmoidA_neg, [sigmoidB_pos]);
%     else
%         f_plot_prob_colormap(strcat(sample_name, '_neg'), ...
%             sigmoidA, sigmoidB_neg);
%         f_plot_prob_colormap(strcat(sample_name, '_pos'), ...
%             sigmoidA, sigmoidB_pos);
%         f_plot_logratio_hists(sample_name, ...
%             logratios_3day_neg, ...
%             logratios_3day_pos);
%         f_plot_logratio_hists_essub(sample_name, ...
%             logratios_3day_neg, ...
%             logratios_3day_pos, ...
%             logratios_es);
%     end
    
    % write out the results for this one
    text_output = strcat( ...
        num2str(ii), ',', ...
        gene, ',', ...
        sample_name, ',', ...
        num2str(sigmoidB_neg), ',', ...
        num2str(sigmoidB_pos), ',', ...
        num2str(frac_me_neg), ',', ...
        num2str(frac_me_pos), ',', ...
        num2str(kl_3day), ',', ...
        num2str(kl_sign), ',', ...
        num2str(exp_frac_me_from_pOS), '\n');
    fprintf(fileID, text_output);
    
end
fclose(fileID);

k = keys(gene_shifts);
for ii=1:length(gene_shifts)
    f_plot_sigmoid_shift(strcat(k{ii}, '_all'), sigmoidA, ...
            0, gene_shifts(k{ii}));
end

k = keys(gene2idx);
bcenters = -7:0.4:10;
for ii=1:length(k)
    idxlist = gene2idx(k{ii});
    fig = figure();
    prob_neg = [];
    prob_pos = [];
    for jj=1:length(idxlist)
        mean_neg = mean(osr_distributions_neg{idxlist(jj)});
        std_neg  = std(osr_distributions_neg{idxlist(jj)});
        neg = (osr_distributions_neg{idxlist(jj)} - mean_neg) ./ std_neg;
        pos = (osr_distributions_pos{idxlist(jj)} - mean_neg) ./ std_neg;
        [counts_neg, ~] = hist(neg, bcenters);
        [counts_pos, ~] = hist(pos, bcenters);
        hold on
        prob_neg = [prob_neg; counts_neg ./ sum(counts_neg)];
        prob_pos = [prob_pos; counts_pos ./ sum(counts_pos)];
%         plot(bcenters, counts_neg ./ sum(counts_neg), 'k--')
%         plot(bcenters, counts_pos ./ sum(counts_pos), 'b')
    end
    if length(idxlist) > 1
        mean_neg = mean(prob_neg);
        mean_pos = mean(prob_pos);
        std_neg  = std(prob_neg);
        std_pos  = std(prob_pos);
        patch([bcenters fliplr(bcenters)], ...
            [mean_neg - std_neg fliplr(mean_neg + std_neg)], ...
            [0.9 0.9 0.9], ...
            'EdgeColor', 'none')
        patch([bcenters fliplr(bcenters)], ...
            [mean_pos - std_pos fliplr(mean_pos + std_pos)], ...
            [0.9 0.9 1.0], ...
            'EdgeColor', 'none')
        plot(bcenters, mean_neg, 'k')
        plot(bcenters, mean_pos, 'b')
    else
        plot(bcenters, prob_neg, 'k')
        plot(bcenters, prob_pos, 'b')
    end
    saveas(fig, strcat('../plots/logratio_hists/', k{ii}, ...
           '_all_loghists.pdf'))
    close();
end

end