function decoding_EEG(cfg)

if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 6;end
if ~isfield(cfg, 'nch'); cfg.nch = 100;end
if ~isfield(cfg, 'classifier'); cfg.classifier = @cosmo_classify_lda;end
if ~isfield(cfg, 'var_threshold'); cfg.var_threshold = 0.99;end
if ~isfield(cfg, 'plotting'); cfg.plotting = false;end

%% Define Parameters
s_idx=1;
ft_defaults;

%% Decoding
for s=cfg.subNums % for all subjects

    filepath = fullfile(cfg.outputPath, ['sub-', num2str(s)], 'eeg', ['PEP_WP4_EEG', num2str(s), '_timelock_reref_w', '.mat']);
    load(filepath);

    % get time info
    res.time=timelock.time;

    % convert to cosmo
    ds=cosmo_meeg_dataset(timelock);
    clear timelock

    % add targets
    ds.sa.targets=ds.sa.trialinfo(:,1);

    % add chunks
    ds.sa.chunks=[1:length(ds.sa.targets)]';
    ds.sa.chunks=cosmo_chunkize(ds,cfg.nch);

    fprintf(['Subject ' num2str(s_idx) ' of ' num2str(cfg.n) '.']);

    % conduct decoding analysis at every time point
    for tp=1:max(ds.fa.time)% for all time points

%         display(['Subject ' num2str(s_idx) ' of ' num2str(cfg.n) '. Time point ' num2str(tp) ' of ' num2str(max(ds.fa.time)) '.']);

        ds_tp=cosmo_slice(ds,ismember(ds.fa.time,tp),2);

        % do pca
        n_feat=length(unique(ds_tp.fa.chan));
        [coeff,x,LATENT,~,x_exp,mu]=pca(ds_tp.samples);
        for ccx=1:length(x_exp)
            if sum(x_exp(1:ccx))>=cfg.var_threshold*100
                n_feat=ccx;
                break
            end
        end
        ds_tp.samples=x(:,1:n_feat);
        ds_tp.fa.chan=1:n_feat;
        ds_tp.fa.time=repmat(tp,1,n_feat);
        ds_tp.a.fdim.values{1}=1:n_feat;

        ds_class=ds_tp;

        % assign targets (i.e. scene trigger codes) to categories
        % 1-50 bathrooms
        % 51-100 kitchen

        bathroom_logical=ismember(ds.sa.targets, 1:50);
        kitchen_logical=ismember(ds.sa.targets, 51:100);

        ds_class.sa.targets(bathroom_logical,1)=1;
        ds_class.sa.targets(kitchen_logical,1)=2;

        % decoding settings
        args.classifier=cfg.classifier;
        args.partitions=cosmo_nchoosek_partitioner(ds_class,1);

        % run decoding
        acc=cosmo_crossvalidation_measure(ds_class,args);
        res.dec_acc(s_idx,tp)=acc.samples;
        res.order(s_idx) = s;

    end % time points
    s_idx=s_idx+1;
end % subjects

% smooth the decoding curve with a rolling average
res.dec_acc=smoothdata(res.dec_acc,2,'movmean',cfg.smoothing_window);

file_name=['PEP_WP4_EEG_decoding_accuracy',num2str(cfg.smoothing_window),'_TP_RA_reref_w.mat'];

% save decoding data to file
if ~exist(cfg.outputPath, 'dir')
    mkdir(cfg.outputPath);
end

cd(cfg.outputPath);
save(file_name,'res');

%% plot decoding results
if cfg.plotting

    % graphic parameters
    set(0, 'DefaultTextFontSize', cfg.FontSize);
    set(0, 'DefaultAxesFontSize', cfg.FontSize);
    set(0, 'DefaultTextFontName', cfg.FontName);
    set(0, 'DefaultAxesFontName', cfg.FontName);

    dec_acc_mean=squeeze(mean(res.dec_acc,1));

    fig=figure();
    plot(res.time, dec_acc_mean(:), 'LineWidth', 2);
    yline(1/2,'k--', 'LineWidth', 2, 'Color', 'r');
    text(0.5, 0.5, 'Chance', 'Color', 'r', ...
        'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'left');
    yline(max(dec_acc_mean), 'k--', ['max ', num2str(max(dec_acc_mean))])
    xline(0,'k', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Decoding accuracy');
    title('Decoding accuracy');
    set(gcf, 'color', [1 1 1]); % white background
    set(gca, 'box', 'off');

    if cfg.saving
        save_plot(fig, 'decoding_accuracy_reref_w', cfg.figPath);
    end
end
end