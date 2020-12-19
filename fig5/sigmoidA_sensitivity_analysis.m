close all;
clear all;

rng(1)

sigmoidAs = linspace(1,10,91);

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
    
    % output file location
    floc = strcat('results/results_', num2str(sigmoidA), '.csv');
    if isfile(floc)
        continue
    end
    fileID = fopen(floc, 'w');
    
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
        frac_me_neg = f_get_frac_me(data_endpt_neg, coords_me, ...
            coords_ecto);
        frac_me_pos = f_get_frac_me(data_endpt_pos, coords_me, ...
            coords_ecto);
        
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
        end
        
        % calculate log ratios
        logratios_3day_neg = log2(data_3day_neg.Red ...
            ./ data_3day_neg.Green);
        logratios_3day_pos = log2(data_3day_pos.Red ...
            ./ data_3day_pos.Green);
        
        % infer where the sigmoid is
        sigmoidB_neg = ...
            f_infer_sigmoid_loc(logratios_3day_neg, frac_me_neg, sigmoidA);
        sigmoidB_pos = ...
            f_infer_sigmoid_loc(logratios_3day_pos, frac_me_pos, sigmoidA);
        if ~isKey(gene_shifts, gene)
            gene_shifts(gene) = [];
        end
        gene_shifts(gene) = [gene_shifts(gene), ...
            sigmoidB_pos - sigmoidB_neg];
        
        % store the OSR distributions
        osr_distributions_neg{ii} = logratios_3day_neg;
        osr_distributions_pos{ii} = logratios_3day_pos;
        
        % calculate the KL divergence between the cfp+ and cfp- 
        % distributions at 3 days
        kl_3day = f_kl_div(logratios_3day_neg, logratios_3day_pos, ...
            nbins, pseudocount);
        kl_sign = NaN;
        if mean(logratios_3day_pos) > mean(logratios_3day_neg)
            kl_sign = 1;
        else
            kl_sign = -1;
        end
        
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
end

floc = strcat('results/results_8.csv');
        fileID = fopen(floc);
        m = textscan(fileID, ...
            '%d %s %s %f %f %f %f %f %d %f', ...
            'Delimiter', ',', ...
            'HeaderLines', 1);
        fclose(fileID);
genes = m{2};


for gg=1:length(genes)
    clf;
    saveit = true;
    gene = genes{gg};
    means = [];
    stdevs = [];
    sigs = [];
    for aa=1:length(sigmoidAs)
        sigmoidA = sigmoidAs(aa);
        floc = strcat('results/results_', num2str(sigmoidA), '.csv');
        fileID = fopen(floc);
        m = textscan(fileID, ...
            '%d %s %s %f %f %f %f %f %d %f', ...
            'Delimiter', ',', ...
            'HeaderLines', 1);
        fclose(fileID);
        
        sigmoidB_neg_cfp = m{4}(strcmp(m{2}, 'CFP'));
        sigmoidB_pos_cfp = m{5}(strcmp(m{2}, 'CFP'));
        cfp_shifts = sigmoidB_pos_cfp - sigmoidB_neg_cfp;
        cfp_shift = mean(cfp_shifts);
        
        sigmoidB_neg = m{4}(strcmp(m{2}, gene));
        if length(sigmoidB_neg) < 3
            saveit = false;
            break
        end
        sigmoidB_pos = m{5}(strcmp(m{2}, gene));
        [~, pval] = ttest2(sigmoidB_pos-sigmoidB_neg, cfp_shifts);
        sigs = [sigs, pval < 0.05];
        sigmoid_shift = sigmoidB_pos - sigmoidB_neg - cfp_shift;
        shift_mean = mean(sigmoid_shift);
        means = [means shift_mean];
        shift_stdev = std(sigmoid_shift);
        stdevs = [stdevs shift_stdev];
    end
    if saveit
        figure(1)
        hold on
        ymin = min([0, means + stdevs, means - stdevs]);
        ymax = max([0, means + stdevs, means - stdevs]);
        patch([sigmoidAs, fliplr(sigmoidAs)], ...
            [means - stdevs, fliplr(means + stdevs)], 0, ...
            'FaceColor', '#CCCCCC', ...
            'EdgeColor', '#FFFFFF');
        plot(sigmoidAs, means, 'k-')
        scatter(sigmoidAs, means, 20, int8(sigs), 'filled')
        r = ymax - ymin;
        ylim([ymin - 0.05 * r, ymax + 0.05 * r])
        caxis([0,1])
        plot([sigmoidAs(1), sigmoidAs(end)], [0, 0], 'k--')
        saveas(gcf, strcat('plots/sigmoidA_sensitivity/', gene, '.png'))
    end
end
close();

