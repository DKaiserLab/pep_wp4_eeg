function res = doPermutations(matin2D, rval, cfg)

% get default values
if ~isfield(cfg, 'permutation_type'); cfg.permutation_type = 'row_col_shuffle_ref';end
if ~isfield(cfg, 'n_permutations'); cfg.n_permutations = 10000;end
if ~isfield(cfg, 'partial_cor'); cfg.partial_cor = false;end
if ~isfield(cfg, 'random_seqs')
    % generate permutated subjects list
    rng(1) % ensure reproducible outcome
    cfg.random_seqs = cell(width(matin2D), cfg.n_permutations);
    for i = 1:cfg.n_permutations
        if ismember(cfg.permutation_type, {'row_col_shuffle_ref',...
                'row_col_shuffle_pred', 'row_col_shuffle_pred_all', 'row_col_shuffle_pred_plus_reoder'})
            % get random sequence of rows and columns
            for ii = 1:width(matin2D)
                cfg.random_seqs{ii,i} = randperm(cfg.n);
            end
        elseif strcmp(cfg.permutation_type, 'sign_flip_ref')
            % get random rows and columns to flip sign
            cfg.random_seqs{1,i} = randperm(cfg.n, randi([1, cfg.n]));
        end
    end
end

% make random permutations
for iRDM = 1:width(matin2D)
    permutation_RDMs(iRDM).RDM = squareform(matin2D(:, iRDM));
    permutation_RDMs(iRDM).name = ['RDM', num2str(iRDM)];
end
orginal_RDMs = permutation_RDMs;
res.perm_r_mat = zeros(width(matin2D)-1, cfg.n_permutations);

for perm = 1:cfg.n_permutations
    % randomize RDM
    if strcmp(cfg.permutation_type, 'row_col_shuffle_ref')
        % shuffle the order of rows and columns
        ref_RDM = orginal_RDMs(1);
        ref_RDM.RDM = ref_RDM.RDM(cfg.random_seqs{1,perm}, cfg.random_seqs{1,perm});
        % replace reference RDM by permutated RDM in RDMs struct
        permutation_RDMs(1) = ref_RDM;
    elseif strcmp(cfg.permutation_type, 'sign_flip_ref')
        ref_RDM_rand = ref_RDM;
        % flip sign of row
        ref_RDM_rand.RDM(cfg.random_seqs{1,perm},:) = -(ref_RDM.RDM(cfg.random_seqs{1,perm},:));
        % flip sign of columns
        ref_RDM_rand.RDM(:,cfg.random_seqs{1,perm}) = -(ref_RDM.RDM(:,cfg.random_seqs{1,perm}));
        % replace reference RDM by permutated RDM in RDMs struct
        permutation_RDMs(1) = ref_RDM_rand;
    end

    % run  correlation
    if cfg.partial_cor
        [~, r_mat, ~, ~] = partial_cor_RDM(cfg, permutation_RDMs);
    else
        [~, r_mat, ~] = cor_RDM(permutation_RDMs, cfg);
    end
    res.perm_r_mat(1:end, perm) = r_mat(2:end, 1);


    if mod(perm/cfg.n_permutations, 0.1) == 0
        disp([num2str((perm/cfg.n_permutations)*100), '% of permutations of is done'])
    end
end
% get confidence interval and p value
% get p value from random permutation (one-sided test aginst
% permutation distribution)
for row = height(res.perm_r_mat)
    res.p_value(row) = sum(res.perm_r_mat(row,:) >= rval(row)) / cfg.n_permutations;

    % get p values for categories
    res.ci_upper(row) = rval - prctile(res.perm_r_mat(row,:), 5);
    res.ci_lower(row) = rval - prctile(res.perm_r_mat(row,:), 95);
end
end