function stats = stats_permutationtest(cfg, d)

% evaluate input
if ~isfield(cfg, 'RDM_to_partial_out'); cfg.RDM_to_partial_out = {'typical_late', 'control_late'}; end
if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'pearson'; end
if ~isfield(cfg, 'add_legend'); cfg.add_legend = true; end
if ~isfield(cfg, 'show_single_cate'); cfg.show_single_cate = false; end
if ~isfield(cfg, 'partial_cor'); cfg.partial_cor = true; end
if ~isfield(cfg, 'save_name'); cfg.save_name = 'compare_roi_RDMs_to_predictor_RDMs'; end
if ~isfield(cfg, 'xaxis_labels'); cfg.xaxis_labels = true; end
if ~isfield(cfg, 'dnns'); cfg.dnns = {cfg.dnn}; end
if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 6; end
if ~isfield(cfg, 'n_permutations'); cfg.n_permutations = 5000; end
if ~isfield(cfg, 'use_maxT'); cfg.use_maxT = false; end
if ~isfield(cfg, 'alpha'); cfg.alpha = 0.05; end

cfg.plot_rdm = false;
cfg.ISC_type = 'pairRep';
partial = cfg.partial_cor;

n_tp = length(d.(cfg.ISC_type).included_time);
pos_tp = d.pairRep.included_time(d.pairRep.all_time(d.pairRep.included_time)>=0);
n_pos_tp = length(pos_tp);
n_predictors = length(cfg.RDM_to_partial_out);
n_sub = cfg.n;
predictors = cfg.RDM_to_partial_out;
idx = 0;

% init
pvalMat = struct;
pvalMat.twosided = struct;
pvalMat.onesided = struct;
percBounds = struct;
sigMat = struct;

% rng(1);

% loop tp
for iTp = 1:n_tp
    act_tp = d.pairRep.all_time(d.pairRep.included_time(iTp));

    if act_tp < 0
        % add nan
        continue;
    end

    idx = idx + 1;

    % init data storage
    perm_corrs_allCats = cell(length(cfg.categories),1);
    r_obs_allCats = cell(length(cfg.categories),1);

    % loop cat
    for cate_num = 1:length(cfg.categories)

        category = char(cfg.categories{cate_num});

        % init RDM
        RDMs = struct;
        RDMs(1).name = d.ISC.([category,'_RDM']).(cfg.ISC_type).name;
        RDMs(1).color = d.ISC.([category,'_RDM']).(cfg.ISC_type).color;
        RDMs(1).RDM = d.ISC.([category,'_RDM']).(cfg.ISC_type)(1).RDM;

        labels = {RDMs.name};
        [RDMs, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, category);

        RDM_orig = RDMs(1).RDM; % save original once per timepoint

        %
        r_obs = d.resMat.partial_cor.(category)(:, iTp);
        r_obs_allCats{cate_num} = r_obs;

        % init data storage
        perm_corrs = nan(n_predictors, cfg.n_permutations);


        if isempty(gcp('nocreate'))
            parpool(8);
        end

        % create permutations and compute correlations
        %parfor
        parfor p = 1:cfg.n_permutations
            s = RandStream('mrg32k3a','Seed', 1 + p);

            perm_idx = randperm(s, n_sub);

            RDMs_local = RDMs;
            RDMs_local(1).RDM = RDM_orig(perm_idx, perm_idx);

            %             perm_RDM = RDM_orig(perm_idx, perm_idx); % permute original, shuffle rows and columns
            %             RDMs(1).RDM = perm_RDM;

            if partial
                [~, rMat_perm, ~] = partial_cor_RDM(cfg, RDMs_local);
            else
                [~, rMat_perm, ~] = cor_RDM(RDMs_local, cfg);
            end

            perm_corrs(:, p) = rMat_perm(2:end, 1);
        end

        nUnique = size(unique(perm_corrs','rows'),1);
        fprintf('Permutations: %d requested, %d unique\n', cfg.n_permutations, nUnique);
        fprintf('Varianz Perm corr (mean across preds): %g\n', mean(var(perm_corrs,0,2)));


        perm_corrs_allCats{cate_num} = perm_corrs;

        for v = 1:n_predictors
            if act_tp == d.pairRep.all_time(pos_tp(1)) && v == 1
                pvalMat.twosided.(category) = nan(n_predictors, n_pos_tp);
                pvalMat.onesided.(category) = nan(n_predictors, n_pos_tp);
                permDistributions.(category) = cell(n_predictors, n_pos_tp);
                percBounds.twosided.(category) = struct;
                percBounds.twosided.(category).upper = nan(n_predictors, n_pos_tp);
                percBounds.twosided.(category).lower = nan(n_predictors, n_pos_tp);
                percBounds.onesided.(category) = struct;
                percBounds.onesided.(category).upper = nan(n_predictors, n_pos_tp);
            end
            % save permutations
            permDistributions.(category){v, idx} = perm_corrs(v, :);
            % compute and save p-values
            pvalMat.twosided.(category)(v, idx) = mean(abs(perm_corrs(v,:)) >= abs(r_obs(v))); % two-sided
            pvalMat.onesided.(category)(v, idx) = mean(perm_corrs(v,:) >= r_obs(v)); % one-sided;
            % compute and save percentiles
            percBounds.twosided.(category).upper(v, idx) = prctile(perm_corrs(v,:), 100-(cfg.alpha*100/2));
            percBounds.twosided.(category).lower(v, idx) = prctile(perm_corrs(v,:), 0+(cfg.alpha*100/2));
            percBounds.onesided.(category).upper(v, idx) = prctile(perm_corrs(v,:), 100-cfg.alpha*100);

        end

    end % category

    r_obs_mean = mean(cat(2, r_obs_allCats{:}),2);

    perm_corrs_mean = nan(n_predictors, cfg.n_permutations);
    for p = 1:cfg.n_permutations
        vals_per_cat = cellfun(@(x) x(:,p), perm_corrs_allCats, 'UniformOutput', false);
        perm_corrs_mean(:,p) = mean(cat(2, vals_per_cat{:}), 2);
    end

    for v = 1:n_predictors

        if act_tp == d.pairRep.all_time(pos_tp(1)) && v == 1
            pvalMat.twosided.all = nan(n_predictors, n_pos_tp);
            pvalMat.onesided.all = nan(n_predictors, n_pos_tp);
            permDistributions.all = cell(n_predictors, n_pos_tp);
            percBounds.twosided.all = struct;
            percBounds.twosided.all.upper = nan(n_predictors, n_pos_tp);
            percBounds.twosided.all.lower = nan(n_predictors, n_pos_tp);
            percBounds.onesided.all = struct;
            percBounds.onesided.all.upper = nan(n_predictors, n_pos_tp);
        end

        % save permutations
        permDistributions.all{v, idx} = perm_corrs_mean(v,:);

        % compute and save p-values
        pvalMat.twosided.all(v, idx) = mean(abs(perm_corrs_mean(v,:)) >= abs(r_obs_mean(v))); % two-sided;
        pvalMat.onesided.all(v, idx) = mean(perm_corrs_mean(v,:) >= r_obs_mean(v)); % one-sided

        % compute and save percentiles
        percBounds.twosided.all.upper(v, idx) = prctile(perm_corrs_mean(v,:), 100-(cfg.alpha*100/2));
        percBounds.twosided.all.lower(v, idx) = prctile(perm_corrs_mean(v,:), 0+(cfg.alpha*100/2));
        percBounds.onesided.all.upper(v, idx) = prctile(perm_corrs_mean(v,:), 100-cfg.alpha*100);
    end

end % tp

if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end


pvalMat.onesided.FDR = struct;
pvalMat.onesided.FDR.all = nan(n_predictors, n_pos_tp);
sigMat.FDR = struct;
sigMat.FDR.all = false(n_predictors, n_pos_tp);

% FDR, separate for each modell
for i=1:n_predictors
    pvals = pvalMat.onesided.all(i,:);
    [is_significant, ~, ~, adj_p] = fdr_bh(pvals, cfg.alpha, 'pdep', 'yes');
    sigMat.FDR.all(i,:) = is_significant;
    pvalMat.onesided.FDR.all(i,:) = adj_p;
end

% FDR, separate for each modell and category
for c = 1:length(cfg.categories)
    cat_name = cfg.categories{c};
    pvalMat.onesided.FDR.(cat_name) = nan(n_predictors, n_pos_tp);
    sigMat.FDR.(cat_name) = false(n_predictors, n_pos_tp);
    for i=1:n_predictors
%         pred = cfg.RDM_to_partial_out{i};
        pvals = pvalMat.onesided.(cat_name)(i,:);
        [is_significant, ~, ~, adj_p] = fdr_bh(pvals, cfg.alpha, 'pdep', 'yes');
        sigMat.FDR.(cat_name)(i,:) = is_significant;
        pvalMat.onesided.FDR.(cat_name)(i,:) = adj_p;
    end
end

stats = struct;
stats.pvalMat = pvalMat;
stats.sigMat = sigMat;
stats.percBounds = percBounds;
stats.predictors = predictors;
stats.testedTime = pos_tp;

end