function [stats, r_vals] = compare_to_drawings_tp(cfg, d)

% evaluate input
if ~isfield(cfg, 'RDM_to_partial_out'); cfg.RDM_to_partial_out = {'typical_late', 'control_late'}; end
if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'pearson'; end
if ~isfield(cfg, 'partial_cor'); cfg.partial_cor = true; end
if ~isfield(cfg, 'n_permutations'); cfg.n_permutations = 5000; end
if ~isfield(cfg, 'cluster'); cfg.cluster = true; end
if ~isfield(cfg, 'ISC_type'); cfg.ISC_type = 'pairRep'; end
if ~isfield(cfg, 'FDR'); cfg.FDR = true; end
correlation_type = cfg.correlation_type;
plottingPreds = 1;

cfg.plot_rdm = false;
partial = cfg.partial_cor;


% loop over frequency bands
for frq = 1:length(cfg.frequencies)
    frqBand = cfg.frequencies{frq};

    %% init
    tps = d.(cfg.ISC_type).(frqBand).included_time;
    n_tp = length(tps);
    tp_vals = d.(cfg.ISC_type).(frqBand).all_time;
    pos_tp = tps(tp_vals(tps) > 0);
    n_neg_tp = n_tp - length(pos_tp);
    n_predictors = length(cfg.RDM_to_partial_out);
    n_sub = cfg.n;
    all_preds = nan(nchoosek(cfg.n,2), length(cfg.RDM_to_partial_out));

    r_obs = nan(n_tp, plottingPreds, length(cfg.categories));
    perm_corrs = nan(n_tp, cfg.n_permutations , plottingPreds, length(cfg.categories));


    %% start test
    % create permutation matrix
    rng(1);
    perm_idx_mat = zeros(n_sub, cfg.n_permutations);
    for p = 1:cfg.n_permutations
        perm_idx_mat(:,p) = randperm(n_sub);
    end

    if isempty(gcp('nocreate'))
        parpool(8);
    end

    % loop cat
    for c = 1:length(cfg.categories)
        category = char(cfg.categories{c});

        % get predictors and regress out control predictors
        RDMs = d.DNN.(cfg.dnn).control.(category).subject_mean(1);
        labels = {RDMs.name};
        [RDMs, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, category);
        for iPred = 1:n_predictors
            RDMs(iPred+1).RDM(eye(cfg.n) == 1) = 0;
            all_preds(:, iPred) = squareform(RDMs(iPred+1).RDM);
        end
        X = [all_preds(:,2:end), ones(size(all_preds,1),1)];
        b = X \ all_preds(:,1);
        uniquePred = all_preds(:,1) - X*b;

        %% loop tp
        dispstat('','init')
        for iTp = 1:n_tp
            act_tp = tp_vals(tps(iTp));

            % progress report
            if mod(iTp, 10) == 0
                dispstat(['TP ', num2str(iTp), ' out of ', num2str(n_tp),...
                    ' - ', category, ' - ', frqBand])
            end

            % init RDM
            eegRDM = d.ISC.([category,'_RDM']).(frqBand).(cfg.ISC_type)(iTp).RDM;
            eegRDM(eye(size(eegRDM)) == 1) = 0;
            r_obs(iTp, 1, c) = corr(squareform(eegRDM)', uniquePred, 'row', 'pairwise', 'type', correlation_type);


            % do permutationtest only on time points after stimulus onset
            if act_tp > 0

                % create permutations and compute correlations
                allPerms = nan(cfg.n_permutations, length(squareform(eegRDM)));
                rng(1)
                parfor p = 1:cfg.n_permutations

                    % shuffle neural RDM
                    RDM_shuffled = eegRDM(perm_idx_mat(:,p), perm_idx_mat(:,p));
                    allPerms(p, :) = squareform(RDM_shuffled);
                end
                perm_corrs(iTp, :, 1, c) = corr(allPerms', uniquePred, 'row', 'pairwise', 'type', correlation_type);
            end

        end % tp
    end % category

    % take mean of categories 
    r_obs_mean = mean(r_obs, 3);
    perm_corrs_mean = mean(perm_corrs, 4);

    %% stats

    % p values
    stats.(frqBand).p_vals = nan(size(r_obs_mean));
    stats.(frqBand).ci = nan(length(r_obs_mean), 2);
    for iTp = 1:length(r_obs_mean)

        % only for postive time points
        if tp_vals(tps(iTp)) > 0
            stats.(frqBand).p_vals(iTp) = sum(perm_corrs_mean(iTp, :) >= r_obs_mean(iTp))/cfg.n_permutations; % p value (one sided)
            stats.(frqBand).ci(iTp, 1) = prctile(perm_corrs_mean(iTp, :), 95); % upper ci
            stats.(frqBand).ci(iTp, 2) = prctile(perm_corrs_mean(iTp, :), 5); % lower ci
        end
    end

    % fdr correction
    stats.(frqBand).fdr_p = nan(size(r_obs_mean));
    noNan = ~isnan(stats.(frqBand).p_vals);
    [~, ~, ~, stats.(frqBand).fdr_p(noNan)] = ...
        fdr_bh(stats.(frqBand).p_vals(noNan));

    % cluster
    stats.(frqBand).cluster.all = stats_cluster(r_obs_mean, ...
        squeeze(perm_corrs_mean), ...
        cfg, d, frqBand);

    stats.(frqBand).cluster_p = ones(size(r_obs_mean));
    if ~isempty(stats.(frqBand).cluster.all.pvals_cluster)
        cluster_points = stats.(frqBand).cluster.all.obs_clusters.PixelIdxList{:};
        stats.(frqBand).cluster_p(cluster_points) = min(stats.(frqBand).cluster.all.pvals_cluster);
    end

    % write to data struct
    r_vals.(frqBand).r_obs_mean = r_obs_mean;

end

end