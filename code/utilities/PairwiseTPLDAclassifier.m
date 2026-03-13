function [res, meanAcc, cfg] = PairwiseTPLDAclassifier(cfg)

% evaluate input
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'pca'); cfg.pca = false; end
if ~isfield(cfg, 'zscore'); cfg.zscore = true; end
if ~isfield(cfg, 'frequencies'); cfg.frequencies = {'full'}; end
if ~isfield(cfg, 'decoding_start'); cfg.decoding_start = -0.2; end
if ~isfield(cfg, 'decoding_end'); cfg.decoding_end = 0.5; end
if ~isfield(cfg, 'makeBetweenComparison'); cfg.makeBetweenComparison = false; end
makeBetweenComparison = cfg.makeBetweenComparison;
if ~isfield(cfg, 'nTrials2average'); cfg.nTrials2average = 1; end
if cfg.nTrials2average ~= 1 && mod(cfg.nTrials2average, 2) ~= 0
    error('nTrials2average must be 1 or a even number')
end
if ~isfield(cfg, 'nBlocks'); cfg.nBlocks = 20; end


% Define classifiers
classifier = @cosmo_classify_lda;

% Define frequency bands
cfg.frqBands.alpha = [8, 12];
cfg.frqBands.beta = [13, 30];
cfg.frqBands.gamma = [31, 70];

% loop over frequency bands
for frq = 1:length(cfg.frequencies)
    frqBand = cfg.frequencies{frq};

    % get filename
    freqTag = [];
    if ~strcmp(frqBand, 'full')
        freqTag = ['_', frqBand];
    end

    avgTag = [];
    if cfg.nTrials2average > 1
        avgTag = '_avg';
    end

    pcaTag = [];
    if cfg.pca
        pcaTag = '_pca';
    end

    % loop through subjects
    for iSub = 1:length(cfg.subNums)

        subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
        subID2 = strrep(subID, '-', '');

        % make file name
        fileName = fullfile(cfg.outputPath, subID, 'eeg',...
            ['PEP_WP4_EEG_', num2str(cfg.subNums(iSub)), '_pairwise_decoding'...
            freqTag, avgTag, pcaTag, ...
            '.mat']);

        % load data if exist
        if exist(fileName, 'file')
            disp(['Loading pairwise decoding for subject ',  subID]);
            load(fileName)
        else


            % progress report
            disp(['Starting pairwise decoding for subject ',  subID]);

            % get preprocessed data
            if strcmp(frqBand, 'full')
                filepath = fullfile(cfg.outputPath, ['sub-', num2str(cfg.subNums(iSub))], 'eeg',...
                    ['PEP_WP4_EEG', num2str(cfg.subNums(iSub)), '_timelock_reref_s2', '.mat']);
            else
                currentBand = cfg.frqBands.(frqBand);
                filepath = fullfile(cfg.outputPath, ['sub-', num2str(cfg.subNums(iSub))], 'eeg',...
                    ['PEP_WP4_EEG', num2str(cfg.subNums(iSub)), '_tfa', '.mat']);
            end
            load(filepath);

            %convert to cosmo
            if ~strcmp(frqBand, 'full')
                timelock = tf; clear tf
            end
            ds=cosmo_meeg_dataset(timelock);

            % freq band for decoding
            if ~strcmp(frqBand, 'full')

                % filter current frequency band
                filt_band = (ds.fa.freq >= currentBand(1)) & (ds.fa.freq <= currentBand(2));
                ds = cosmo_slice(ds, filt_band, 2);
            end

            % check if data set is complete
            if height(ds.samples) == cfg.nTrials*cfg.nBlocks

                % define chunks
                nch=cfg.nBlocks/cfg.nTrials2average; % is based on many repetitions we have and over how many trials we average
                %nch = 20;
                ds.sa.chunks=(1:length(ds.sa.trialinfo))';
                ds.sa.targets=ds.sa.trialinfo;
                ds.sa.chunks=cosmo_chunkize(ds, nch);

                % averaging trials
                ds = cosmo_average_samples(ds, 'targets', ds.sa.targets, 'chunks', ds.sa.chunks);

            else
                warning(['Subject: ', num2str(cfg.subNums(iSub)), ' has incomplete dataset'])
                ds = makeChunkForIncompleteDS(ds, cfg, iSub);
            end

            % time range for decoding
            timepoints = find(timelock.time >= cfg.decoding_start &...
                timelock.time <= cfg.decoding_end);

            %get time info
            all_time=timelock.time;
            included_time=timepoints;
            clear timelock %new

            %% Pairwise decoding
            % Initialize RDM
            rdm = zeros(cfg.nTrials, cfg.nTrials, length(timepoints));
            mean_accuracy = zeros(1, length(timepoints));
            disp('')

            % remove nans
            ds = cosmo_remove_useless_data(ds);

            if cfg.zscore
                ds = cosmo_normalize(ds, 'zscore');
            end

            dispstat('','init')
            for tp=1:length(timepoints)%1:max(ds.fa.time)

                if mod(tp, 10) == 0
                    dispstat(['TP ', num2str(tp), ' out of ', num2str(length(timepoints))])
                end

                % select time point
                tp_idx = timepoints(tp);
                ds_tp=cosmo_slice(ds, ismember(ds.fa.time,tp_idx),2);

                if cfg.pca
                    % do PCA
                    [ds_pca, ~] = cosmo_map_pca(ds_tp, 'pca_explained_ratio', 0.9, 'max_feature_count', 5000);
                    ds_tp.samples = ds_pca.samples;
                    ds_tp.fa = ds_pca.fa;
                    ds_tp.a = ds_pca.a;
                end

               
                if isempty(gcp('nocreate'))
                    parpool(10);
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

            % Save decoding resultst
            save(fileName, "rdm", "mean_accuracy", "all_time", "included_time")

        end

        % Write RDM to struct
        res.(subID2).(frqBand).rdm = rdm;
        res.(subID2).(frqBand).mean_accuracy = mean_accuracy;


    end % sub

    for cate_num = 1:numel(cfg.categories)
        category = char(cfg.categories{cate_num});

        % preallocate
        meanAcc_temp = nan(length(included_time), cfg.n);

        for tp=1:length(included_time)
            for iSub = 1:length(cfg.subNums)
                subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
                subID2 = strrep(subID, '-', '');

                % make a matrix with vectorized RDMs
                rdm = squeeze(res.(subID2).(frqBand).rdm(:,:,tp));
                rdm(eye(size(rdm)) == 1) = 0;

                % filter for category
                if strcmp(category, cfg.categories{1})
                    rdm = rdm(1:cfg.nTrials/2, 1:cfg.nTrials/2);
                else
                    rdm = rdm(cfg.nTrials/2 + 1:end, cfg.nTrials/2 + 1:end);
                end

                meanAcc_temp(tp, iSub) = mean(squareform(rdm), 'omitnan');

            end

        end
        % store in struct
        meanAcc.(category).(frqBand) = meanAcc_temp;


    end

    % store info in res
    res.(frqBand).included_time = included_time;
    res.(frqBand).all_time = all_time;
    res.(frqBand).nTrials2average = cfg.nTrials2average;
    res.(frqBand).pca = cfg.pca;

    % save results
    outputFolder = fullfile(cfg.outputPath, 'group_level', 'RDM');
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    save(fullfile(outputFolder, 'newDecoding'), 'res')

end
end