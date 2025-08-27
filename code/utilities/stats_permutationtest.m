function d = stats_permutationtest(cfg, d)

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

n_tp = length(d.(cfg.ISC_type).included_time);
n_predictors = length(cfg.RDM_to_partial_out);

% init
d.pvalMat.avg = nan(n_predictors, n_tp);
d.pvalMat_onesided.avg = nan(n_predictors, n_tp);
d.permDistributions.all = cell(n_predictors, n_tp);

% loop tp
for iTp = 1:n_tp

    % init data storage
    perm_corrs_allCats = cell(length(cfg.categories),1);
    r_obs_allCats = cell(length(cfg.categories),1);

    % loop cat
    for cate_num = 1:length(cfg.categories)

        rng(1)

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
        perm_corrs = zeros(n_predictors, cfg.n_permutations);

        % create permutations and compute correlations
        for p = 1:cfg.n_permutations
            perm_idx = randperm(cfg.n);

            perm_RDM = RDM_orig(perm_idx, perm_idx); % permute original, shuffle rows and columns
            RDMs(1).RDM = perm_RDM;

            if cfg.partial_cor
                [~, rMat_perm, ~] = partial_cor_RDM(cfg, RDMs);
            else
                [~, rMat_perm, ~] = cor_RDM(RDMs, cfg);
            end

            perm_corrs(:, p) = rMat_perm(2:end, 1);
        end

        perm_corrs_allCats{cate_num} = perm_corrs;

        for var = 1:n_predictors
            % save permutations
            d.permDistributions.(category){var, iTp} = perm_corrs(var, :);
            % compute and save p-values
            d.pvalMat.(category)(var, iTp) = mean(abs(perm_corrs(var,:)) >= abs(r_obs(var))); % two-sided
            d.pvalMat_onesided.(category)(var, iTp) = mean(perm_corrs(var,:) >= r_obs(var)); % one-sided;
            % compute and save percentiles
            d.percBounds.(category).upper(var, iTp) = prctile(perm_corrs(var,:), 100-(cfg.alpha*100/2));
            d.percBounds.(category).lower(var, iTp) = prctile(perm_corrs(var,:), 0+(cfg.alpha*100/2));
            d.percBounds_onesided.(category).upper(var, iTp) = prctile(perm_corrs(var,:), 100-cfg.alpha*100);
            
        end
        
    end % category

    r_obs_mean = mean(cat(2, r_obs_allCats{:}),2);

    perm_corrs_mean = zeros(n_predictors, cfg.n_permutations);
    for p = 1:cfg.n_permutations
        vals_per_cat = cellfun(@(x) x(:,p), perm_corrs_allCats, 'UniformOutput', false);
        perm_corrs_mean(:,p) = mean(cat(2, vals_per_cat{:}), 2); 
    end

    for var = 1:n_predictors
        % save permutations
        d.permDistributions.all{var, iTp} = perm_corrs_mean(var,:);
        
         % compute and save p-values
        d.pvalMat.avg(var, iTp) = mean(abs(perm_corrs_mean(var,:)) >= abs(r_obs_mean(var))); % two-sided;
        d.pvalMat_onesided.avg(var, iTp) = mean(perm_corrs_mean(var,:) >= r_obs_mean(var)); % one-sided
       
        % compute and save percentiles
        d.percBounds.all.upper(var, iTp) = prctile(perm_corrs_mean(var,:), 100-(cfg.alpha*100/2));
        d.percBounds.all.lower(var, iTp) = prctile(perm_corrs_mean(var,:), 0+(cfg.alpha*100/2));
        d.percBounds_onesided.all.upper(var, iTp) = prctile(perm_corrs_mean(var,:), 100-cfg.alpha*100);
    end

end % tp

% FDR, separate for each modell
for i=1:n_predictors
    pred = cfg.RDM_to_partial_out{i};
    pvals = d.pvalMat_onesided.avg(i,:);
    [is_significant, ~, ~, adj_p] = fdr_bh(pvals, cfg.alpha, 'pdep', 'yes');
    d.sigMat_FDR.avg.(pred) = is_significant;
    d.pval_fdr.avg.(pred) = adj_p;
end

% FDR, separate for each modell and category
for c = 1:numel(cfg.categories)
    cat_name = cfg.categories{c};
    for i=1:n_predictors
        pred = cfg.RDM_to_partial_out{i};
        pvals = d.pvalMat_onesided.(cat_name)(i,:);
        [is_significant, ~, ~, adj_p] = fdr_bh(pvals, cfg.alpha, 'pdep', 'yes');
        d.sigMat_FDR.(pred).(cat_name) = is_significant;
        d.pval_fdr.(pred).(cat_name) = adj_p;
    end
end


end