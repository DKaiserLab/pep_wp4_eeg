function stats = stats_permutationtest(cfg, d)

% evaluate input
if ~isfield(cfg, 'RDM_to_partial_out'); cfg.RDM_to_partial_out = {'typical_late', 'control_late'}; end
if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'pearson'; end
if ~isfield(cfg, 'partial_cor'); cfg.partial_cor = true; end
if ~isfield(cfg, 'n_permutations'); cfg.n_permutations = 5000; end
if ~isfield(cfg, 'alpha'); cfg.alpha = 0.05; end
if ~isfield(cfg, 'cluster'); cfg.cluster = true; end
if ~isfield(cfg, 'FDR'); cfg.FDR = true; end

cfg.plot_rdm = false;
cfg.ISC_type = 'pairRep';
partial = cfg.partial_cor;

n_tp = length(d.(cfg.ISC_type).included_time);
pos_tp = d.pairRep.included_time(d.pairRep.all_time(d.pairRep.included_time)>=0);
n_pos_tp = length(pos_tp);
n_predictors = length(cfg.RDM_to_partial_out);
n_sub = cfg.n;
predictors = cfg.RDM_to_partial_out;

%% init
stats = struct;
pvalMat = struct;
pvalMat.twosided = struct;
pvalMat.twosided.all = nan(n_predictors, n_pos_tp);
pvalMat.onesided = struct;
pvalMat.onesided.FDR = struct;
pvalMat.onesided.FDR.all = nan(n_predictors, n_pos_tp);
pvalMat.onesided.all = nan(n_predictors, n_pos_tp);
sigMat = struct;
sigMat.FDR = struct;
sigMat.FDR.all = false(n_predictors, n_pos_tp);

r_obs_mean_all = nan(n_predictors, n_pos_tp); % save mean observed correlations for each predictor and tp
perm_corrs_mean_all = nan(n_predictors, n_pos_tp, cfg.n_permutations); % save mean correlations of permutations for each predictor and tp
r_obs_byCat = cell(length(cfg.categories),1); % save observed correlations for each predictor and tp
perm_corrs_byCat = cell(length(cfg.categories),1); % save correlations of permutations for each predictor and tp

for c = 1:length(cfg.categories)
    category = char(cfg.categories{c});
    r_obs_byCat{c} = nan(n_predictors, n_pos_tp); % save observed correlations for each predictor and tp (and category)
    perm_corrs_byCat{c} = nan(n_predictors, n_pos_tp, cfg.n_permutations); % save correlations of permutations for each predictor and tp (and category)
    pvalMat.twosided.(category) = nan(n_predictors, n_pos_tp);
    pvalMat.onesided.(category) = nan(n_predictors, n_pos_tp);
    pvalMat.onesided.FDR.(category) = nan(n_predictors, n_pos_tp);
    sigMat.FDR.(category) = false(n_predictors, n_pos_tp);
end

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

idx = 0;
% loop tp
for iTp = 1:n_tp
    act_tp = d.pairRep.all_time(d.pairRep.included_time(iTp));

    % do permutationtest only on time points after stimulus onset
    if act_tp < 0
        continue;
    end

    idx = idx + 1;

    % loop cat
    for c = 1:length(cfg.categories)

        category = char(cfg.categories{c});

        % init RDM
        RDMs = struct;
        RDMs(1).name = d.ISC.([category,'_RDM']).(cfg.ISC_type)(iTp).name;
        RDMs(1).color = d.ISC.([category,'_RDM']).(cfg.ISC_type)(iTp).color;
        RDMs(1).RDM = d.ISC.([category,'_RDM']).(cfg.ISC_type)(iTp).RDM;

        labels = {RDMs.name};
        [RDMs, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, category);

        r_obs_byCat{c}(:, idx) = d.resMat.partial_cor.(category)(:, iTp);

        % init data storage
        perm_corrs = nan(n_predictors, cfg.n_permutations);

        % create permutations and compute correlations
        parfor p = 1:cfg.n_permutations
            perm_idx = perm_idx_mat(:,p);

            RDMs_local = RDMs;
            RDMs_local(1).RDM = RDMs(1).RDM(perm_idx, perm_idx);

            if partial
                [~, rMat_perm, ~] = partial_cor_RDM(cfg, RDMs_local);
            else
                [~, rMat_perm, ~] = cor_RDM(RDMs_local, cfg);
            end
            perm_corrs(:, p) = rMat_perm(2:end, 1);
        end
        
        perm_corrs_byCat{c}(:, idx, :) = perm_corrs;
        
        for v = 1:n_predictors
            % compute and save p-values
            pvalMat.twosided.(category)(v, idx) = mean(abs(perm_corrs(v,:)) >= abs(r_obs_byCat{c}(v, idx))); % two-sided
            pvalMat.onesided.(category)(v, idx) = mean(perm_corrs(v,:) >= r_obs_byCat{c}(v, idx)); % one-sided;
        end

    end % category

    vals = cellfun(@(x) x(:,idx), r_obs_byCat, 'UniformOutput', false);
    r_obs_mean_all(:, idx) = mean(cat(2, vals{:}),2);

    for p = 1:cfg.n_permutations
        vals = cellfun(@(x) x(:,idx,p), perm_corrs_byCat, 'UniformOutput', false);
        perm_corrs_mean_all(:, idx, p) = mean(cat(2, vals{:}), 2);
    end


    for v = 1:n_predictors
        % compute and save p-values
        pvalMat.twosided.all(v, idx) = mean(abs(perm_corrs_mean_all(v, idx, :)) >= abs(r_obs_mean_all(v, idx))); % two-sided;
        pvalMat.onesided.all(v, idx) = mean(perm_corrs_mean_all(v, idx, :) >= r_obs_mean_all(v, idx)); % one-sided
    end

end % tp

% if ~isempty(gcp('nocreate'))
%     delete(gcp('nocreate'));
% end

%% compute cluster and/or FDR

% separate for each model
for v = 1:n_predictors
    % FDR
    if cfg.FDR
        pvals = pvalMat.onesided.all(v,:);
        [is_significant, ~, ~, adj_p] = fdr_bh(pvals, cfg.alpha, 'pdep', 'yes');
        sigMat.FDR.all(v,:) = is_significant;
        pvalMat.onesided.FDR.all(v,:) = adj_p;
    end

    % cluster
    if cfg.cluster
        stats.cluster.all{v} = stats_cluster(r_obs_mean_all(v,:), ...
            squeeze(perm_corrs_mean_all(v,:,:)), ...
            cfg, d);
    end

end

% separate for each model and category
for c = 1:length(cfg.categories)
    cat_name = cfg.categories{c};
    
    for v = 1:n_predictors
        % FDR
        if cfg.FDR
            pvals = pvalMat.onesided.(cat_name)(v,:);
            [is_significant, ~, ~, adj_p] = fdr_bh(pvals, cfg.alpha, 'pdep', 'yes');
            sigMat.FDR.(cat_name)(v,:) = is_significant;
            pvalMat.onesided.FDR.(cat_name)(v,:) = adj_p;
        end

        % cluster
        if cfg.cluster
            stats.cluster.(cat_name){v} = stats_cluster(r_obs_byCat{c}(v,:), ...
                squeeze(perm_corrs_byCat{c}(v,:,:)), ...
                cfg, d);
        end
    end
end


% save output
stats.pvalMat = pvalMat;
stats.predictors = predictors;
stats.testedTime = pos_tp;
stats.sigMat = sigMat;


end