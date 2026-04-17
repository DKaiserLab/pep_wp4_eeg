function d = patternISC(d, cfg)

% evaluate input
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'frqFilter'); cfg.frqFilter = false; end
if ~isfield(cfg, 'zscore'); cfg.zscore = true; end
if ~isfield(cfg, 'frequencies'); cfg.frequencies = {'full'}; end
if ~isfield(cfg, 'pISC_start'); cfg.pISC_start = -0.2; end
if ~isfield(cfg, 'pISC_end'); cfg.pISC_end = 0.5; end
if ~isfield(cfg, 'nBlocks'); cfg.nBlocks = 20; end
if ~isfield(cfg, 'regressOutMean'); cfg.regressOutMean = true; end

% Define frequency bands
cfg.frqBands.alpha = [8, 12];
cfg.frqBands.beta = [13, 30];
cfg.frqBands.gamma = [31, 70];

channels = cfg.channels;
nChannels = double(length(channels));

%% loop over frequency bands
for frq = 1:length(cfg.frequencies)
    frqBand = cfg.frequencies{frq};

    % get filename
    freqTag = [];
    if ~strcmp(frqBand, 'full')
        freqTag = ['_', frqBand];
    end

    fltTag = [];
    if cfg.frqFilter
        fltTag = '_filtered';
    end

    regMeanTag = [];
    if cfg.regressOutMean
        regMeanTag = 'regMean';
    end

    % make file name
    folderOut = fullfile(cfg.outputPath, 'group_level', 'pISC');
    fileName = fullfile(folderOut, ...
        ['PEP_WP4_EEG_pISC', freqTag, fltTag, regMeanTag, '.mat']);

    % load data if exist
    if exist(fileName, 'file')
        disp('Loading pISC results');
        load(fileName)


        for tp = 1:size(pISC_mat, 1)

            % Write RDM to struct
            pISC_bat = pISC_mat(tp, 1:cfg.nTrials/2, :);
            pISC_bat_avg = squeeze(mean(pISC_bat));
            d.ISC.bathroom_RDM.(frqBand).pISC(tp).RDM = squareform(pISC_bat_avg);
            d.ISC.bathroom_RDM.(frqBand).pISC(tp).color = [0,0,0];
            d.ISC.bathroom_RDM.(frqBand).pISC(tp).name = num2str(tp);

            pISC_kit = pISC_mat(tp, 1:cfg.nTrials/2, :);
            pISC_kit_avg = squeeze(mean(pISC_kit));
            d.ISC.kitchen_RDM.(frqBand).pISC(tp).RDM = squareform(pISC_kit_avg);
            d.ISC.kitchen_RDM.(frqBand).pISC(tp).color = [0,0,0];
            d.ISC.kitchen_RDM.(frqBand).pISC(tp).name = num2str(tp);

            % store timing
            d.(cfg.ISC_type).(frqBand).included_time = included_time;
            d.(cfg.ISC_type).(frqBand).all_time = all_time;

        end % tp
    else

        disp('Running pISC');

        if isempty(gcp('nocreate'))
            parpool(10);
        end

        %% loop through subjects
        all_ds = {};
        for iSub = 1:length(cfg.subNums)

            subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
            subID2 = strrep(subID, '-', '');


            % progress report
            disp(['Get data for subject ',  subID, ' - ', frqBand]);

            %% get preprocessed data
            if strcmp(frqBand, 'full')
                filepath = fullfile(cfg.outputPath, ['sub-', num2str(cfg.subNums(iSub))], 'eeg',...
                    ['PEP_WP4_EEG', num2str(cfg.subNums(iSub)), '_timelock_reref_s2', fltTag, '.mat']);
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

            % check if data set is complete
            if height(ds.samples) ~= cfg.nTrials*cfg.nBlocks
                warning(['Subject: ', num2str(cfg.subNums(iSub)), ' has incomplete dataset'])
                ds = makeChunkForIncompleteDS(ds, cfg, iSub);
            end

            % averaging all trials with same scene
            ds.sa.targets=ds.sa.trialinfo;
            ds.sa.chunks=ones(length(ds.sa.trialinfo), 1);
            ds = cosmo_average_samples(ds, 'targets', ds.sa.targets, 'chunks', ds.sa.chunks);

            % freq band for decoding
            if ~strcmp(frqBand, 'full')

                % filter current frequency band
                filt_band = (ds.fa.freq >= currentBand(1)) & (ds.fa.freq <= currentBand(2));
                ds = cosmo_slice(ds, filt_band, 2);
            end

            % remove nans
            ds = cosmo_remove_useless_data(ds);

            if cfg.zscore
                ds = cosmo_normalize(ds, 'zscore');
            end

            % store in cell
            all_ds{iSub} = ds;

        end % subs

        % time range for decoding
        timepoints = find(timelock.time >= cfg.pISC_start &...
            timelock.time <= cfg.pISC_end);

        % get time info
        all_time=timelock.time;
        included_time=timepoints;
        clear timelock %new

        % Initialize RDM
        rdm = zeros(cfg.nTrials, cfg.nTrials, length(timepoints));
        mean_accuracy = zeros(1, length(timepoints));
        disp('')

        %% loop through time points
        pISC_mat = nan(length(timepoints), cfg.nTrials, nchoosek(cfg.n, 2));

        dispstat('','init')
        for tp=1:length(timepoints)%1:max(ds.fa.time)

            if mod(tp, 10) == 0
                dispstat(['TP ', num2str(tp), ' out of ', num2str(length(timepoints))])
            end

            % loop through subjects to collect data
            tp_mat = nan(cfg.nTrials,nChannels, cfg.n);
            for iSub = 1:length(cfg.subNums)

                % select time point
                tp_idx = timepoints(tp);
                ds_tp=cosmo_slice(all_ds{iSub}, ismember(all_ds{iSub}.fa.time, tp_idx), 2);
                for c = 1:nChannels
                    c_idx = find(strcmp(channels{c}, ds_tp.a.fdim.values{1}));
                    if ~isempty(c_idx)
                        tp_mat(:, c, iSub) = ds_tp.samples(:, c_idx);
                    end
                end
            end

            % loop through trials
            for tr = 1:cfg.nTrials

                % get ISC
                [~, mat_out, ~] = make_RDM(squeeze(tp_mat(tr, :, :)), cfg);
                mat_out(eye(size(mat_out)) == 1) = 0;
                pISC_mat(tp, tr, :) = squareform(mat_out);

            end % trials

            % Write RDM to struct
            pISC_bat = pISC_mat(tp, 1:cfg.nTrials/2, :);
            pISC_bat_avg = squeeze(mean(pISC_bat));
            d.ISC.bathroom_RDM.(frqBand).pISC(tp).RDM = squareform(pISC_bat_avg);
            d.ISC.bathroom_RDM.(frqBand).pISC(tp).color = [0,0,0];
            d.ISC.bathroom_RDM.(frqBand).pISC(tp).name = num2str(tp);

            pISC_kit = pISC_mat(tp, 1:cfg.nTrials/2, :);
            pISC_kit_avg = squeeze(mean(pISC_kit));
            d.ISC.kitchen_RDM.(frqBand).pISC(tp).RDM = squareform(pISC_kit_avg);
            d.ISC.kitchen_RDM.(frqBand).pISC(tp).color = [0,0,0];
            d.ISC.kitchen_RDM.(frqBand).pISC(tp).name = num2str(tp);

        end % tp

        % Save pISC resultst

        if ~exist(folderOut, 'dir')
            mkdir(folderOut);
        end
        save(fileName, "pISC_mat", "all_time", "included_time")

        % store timing
        d.(cfg.ISC_type).(frqBand).included_time = included_time;
        d.(cfg.ISC_type).(frqBand).all_time = all_time;
    end

end % frg
end
