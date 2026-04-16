function [stats, r_vals] = compare_to_drawings_frq_bands(cfg, d)

% evaluate input

if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'pearson'; end
if ~isfield(cfg, 'partial_cor'); cfg.partial_cor = true; end
if ~isfield(cfg, 'n_permutations'); cfg.n_permutations = 10000; end
if ~isfield(cfg, 'dnn'); cfg.dnn = {'vgg16_imagenet'};end
if ~isfield(cfg, 'dnns'); cfg.dnns = {cfg.dnn};end
if ~isfield(cfg, 'predictor_RDMs'); cfg.predictor_RDMs = {'typical_late', 'control_late', 'photos_late'};end
if ~isfield(cfg, 'RDM_to_partial_out'); cfg.RDM_to_partial_out = cfg.predictor_RDMs; end
if ~isfield(cfg, 'regressOutMean'); cfg.regressOutMean = true; end

correlation_type = cfg.correlation_type;
plottingPreds = 1;
cfg.plot_rdm = false;

% loop over frequency bands
for frq = 1:length(cfg.frequencies)
    frqBand = cfg.frequencies{frq};

    %% init
    n_predictors = length(cfg.RDM_to_partial_out);
    n_sub = cfg.n;
    all_preds = nan(nchoosek(cfg.n,2), length(cfg.RDM_to_partial_out));

    r_obs = nan(plottingPreds, length(cfg.categories));
    perm_corrs = nan(cfg.n_permutations , plottingPreds, length(cfg.categories));

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

        % init RDM
        idx = find(strcmp({d.ISC.kitchen_RDM.(frqBand).blockTfa.name}, 'average'));
        eegRDM = d.ISC.([category,'_RDM']).(frqBand).blockTfa(idx).RDM;
        eegRDM(eye(size(eegRDM)) == 1) = 0;
        r_obs(c) = corr(squareform(eegRDM)', uniquePred, 'row', 'pairwise', 'type', correlation_type);

        % create permutations and compute correlations
        allPerms = nan(cfg.n_permutations, length(squareform(eegRDM)));
        rng(1)
        parfor p = 1:cfg.n_permutations

            % shuffle neural RDM
            RDM_shuffled = eegRDM(perm_idx_mat(:,p), perm_idx_mat(:,p));
            allPerms(p, :) = squareform(RDM_shuffled);
        end
        perm_corrs(:, 1, c) = corr(allPerms', uniquePred, 'row', 'pairwise', 'type', correlation_type);


    end % category

    % take mean of categories
    r_obs_mean = mean(r_obs);
    perm_corrs_mean = mean(perm_corrs, 3);

    %% stats

    % p values
    stats.(frqBand).p_vals = sum(perm_corrs_mean >= r_obs_mean)/cfg.n_permutations; % p value (one sided)
    stats.(frqBand).ci(1) = prctile(perm_corrs_mean, 95); % upper ci
    stats.(frqBand).ci(2) = prctile(perm_corrs_mean, 5); % lower ci
    stats.(frqBand).perm_r = perm_corrs_mean;

    % write to data struct
    r_vals.(frqBand).r_obs_mean = r_obs_mean;

    % display r and p
    disp(newline)
    disp(['r value for ', frqBand, ': ', num2str(r_obs_mean)])
    disp(['p value for ', frqBand, ': ', num2str(stats.(frqBand).p_vals)])

end

end