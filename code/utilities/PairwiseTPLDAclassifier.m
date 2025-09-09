function [res, meanAcc, cfg] = PairwiseTPLDAclassifier(cfg)

% evaluate input
if ~isfield(cfg, 'var_threshold'); cfg.var_threshold = 0.99;end
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'pca'); cfg.pca = false; end
if ~isfield(cfg, 'decoding_start'); cfg.decoding_start = -0.2; end
if ~isfield(cfg, 'decoding_end'); cfg.decoding_end = 0.5; end
if ~isfield(cfg, 'makeBetweenComparison'); cfg.makeBetweenComparison = false; end
makeBetweenComparison = cfg.makeBetweenComparison;
if ~isfield(cfg, 'nTrials2average'); cfg.nTrials2average = 2; end
if cfg.nTrials2average ~= 1 && mod(cfg.nTrials2average, 2) ~= 0
    error('nTrials2average must be 1 or a even number')
end
if ~isfield(cfg, 'nBlocks'); cfg.nBlocks = 20; end

% Define classifiers
classifier = @cosmo_classify_lda;

% get existing pairwise decoding results
fileName = fullfile(cfg.outputPath, 'group_level', 'RDM',...
    'results_RDM_of_pairwise_decoding_reref_avg_pca.mat');
if exist(fileName, 'file')
    load(fileName)
end

% loop through subjects
for iSub = 1:length(cfg.subNums)

    subID = sprintf('sub-%0.3d', cfg.subNums(iSub));

    % check if is part of result structure already
    subID2 = strrep(subID, '-', '');
    if exist('res', 'var')
        if isfield(res, subID2) && cfg.subNums(iSub)~=114
            disp(['RDM for ', subID, ' already exists']);
            continue
        end
    end

    % progress report
    disp(['Starting pairwise decoding for subject ',  subID]);

    filepath = fullfile(cfg.outputPath, ['sub-', num2str(cfg.subNums(iSub))], 'eeg', ['PEP_WP4_EEG', num2str(cfg.subNums(iSub)), '_timelock_reref_w', '.mat']);
    load(filepath);

    %convert to cosmo
    ds=cosmo_meeg_dataset(timelock);

    % check if data set is complete
    if height(ds.samples) == cfg.nTrials*cfg.nBlocks

        % define chunks
        nch=cfg.nBlocks/cfg.nTrials2average; % is based on many repetitions we have and over how many trials we average
        %nch = 20;
        ds.sa.chunks=(1:length(ds.sa.trialinfo))';
        ds.sa.targets=ds.sa.trialinfo;
        ds.sa.chunks=cosmo_chunkize(ds, nch);

    else
        warning(['Subject: ', num2str(cfg.subNums(iSub)), ' has incomplete dataset'])
        ds = makeChunkForIncompleteDS(ds, cfg, iSub);
    end

    % averaging trials
%     ds = cosmo_average_samples(ds);
    ds = cosmo_average_samples(ds, 'targets', ds.sa.targets, 'chunks', ds.sa.chunks);

    % time range for decoding
    timepoints = find(ds.a.fdim.values{2, 1} >= cfg.decoding_start &...
        ds.a.fdim.values{2, 1} <= cfg.decoding_end);

    %get time info
    res.all_time=timelock.time;
    res.included_time=timepoints;
    clear timelock%new

    %% Pairwise decoding
    % Initialize RDM
    rdm = zeros(cfg.nTrials, cfg.nTrials, length(timepoints));
    mean_accuracy = zeros(1, length(timepoints));
    disp('')

    for tp=1:length(timepoints)%1:max(ds.fa.time)

        % select time point
        tp_idx = timepoints(tp);
        ds_tp=cosmo_slice(ds, ismember(ds.fa.time,tp_idx),2);

        if cfg.pca
            % do PCA
            [ds_pca, ~] = cosmo_map_pca(ds_tp,'pca_explained_ratio', cfg.var_threshold);
            ds_tp.samples = ds_pca.samples;
            ds_tp.fa = ds_pca.fa;
            ds_tp.a = ds_pca.a;
        end

        if isempty(gcp('nocreate'))
            parpool(8);
        end
        nTrials = cfg.nTrials;
        parfor stim1 = 1:nTrials
            for stim2 = 1:nTrials
                if ~(stim2 > stim1)
                    continue
                end

                % check whether its a between stimulus comparison
                if ~makeBetweenComparison
                    if (stim2 > nTrials/2 && stim1 <= nTrials/2)
                        rdm(stim2, stim1, tp) = NaN;
                        continue
                    end
                end

                % Subset data for the two stimuli
                ds_stim=cosmo_slice(ds_tp, ds_tp.sa.trialinfo == stim1 | ds_tp.sa.trialinfo == stim2);

                % Rename target
                ds_stim.sa.targets = (ds_stim.sa.trialinfo == stim1) + 1;

                % Define partitions
                partitions = cosmo_nchoosek_partitioner(ds_stim, 1);

                % get predictions for each fold
                [~ ,accuracy] = cosmo_crossvalidate(ds_stim, classifier, partitions);

                % Store the accuracy
                rdm(stim2, stim1, tp) = accuracy;
            end
        end

        % make rdm symetric
        rdm(:, :, tp) = squareform(squareform(rdm(:, :, tp)));

        % take mean
        mean_accuracy(tp) = mean(squareform(rdm(:, :, tp)), 'omitnan');

    end % tp

    % Save the RDM
    res.(subID2).rdm = rdm;
    res.(subID2).mean_accuracy = mean_accuracy;
end % sub

for cate_num = 1:numel(cfg.categories)
    category = char(cfg.categories{cate_num});

    % preallocate
    meanAcc_temp = nan(length(res.included_time), cfg.n);

    for tp=1:length(res.included_time)
        for iSub = 1:length(cfg.subNums)
            subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
            subID2 = strrep(subID, '-', '');

            % make a matrix with vectorized RDMs
            rdm = squeeze(res.(subID2).rdm(:,:,tp));
            rdm(eye(size(rdm)) == 1) = 0;

            % filter for category
            if strcmp(category, cfg.categories{1})
                rdm = rdm(1:cfg.nTrials/2, 1:cfg.nTrials/2);
            else
                rdm = rdm(cfg.nTrials/2 + 1:end, cfg.nTrials/2 + 1:end);
            end

            meanAcc_temp(tp, iSub) = mean(squareform(rdm), 'omitnan');

        end
        meanAcc.(category) = meanAcc_temp;
    end
end



% save results
outputFolder = fullfile(cfg.outputPath, 'group_level', 'RDM');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
save(fileName, 'res')

if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

if cfg.plotting
    PairwiseDecodingPlots(cfg, res);
end

end