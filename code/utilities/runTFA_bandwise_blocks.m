function runTFA_bandwise_blocks(exp_cfg)

if ~isfield(exp_cfg, 'normalize'); exp_cfg.normalize = false; end
if ~isfield(exp_cfg, 'cutResponseWindow'); exp_cfg.cutResponseWindow = true; end

for s = exp_cfg.subNums

    %% -------------------------
    % Paths and filenames
    % -------------------------

    prepro_path = fullfile(exp_cfg.outputPath, ['sub-', num2str(s)], 'eeg');
    tfa_path = prepro_path;
    beh_path = fullfile(exp_cfg.sourcedataPath, ['sub-', num2str(s)], 'beh');

    prepro_filename = ['PEP_WP4_EEG', num2str(s), '_blocks.mat'];
    beh_filename = ['sub-', num2str(s), '_task-main_events.tsv'];
    tfa_filename_all = ['PEP_WP4_EEG', num2str(s), '_tfa_blocks_all.mat'];


    %% -------------------------
    % Skip if exists
    % -------------------------

    if exist(fullfile(tfa_path, tfa_filename_all), 'file')
        fprintf(['File already exists: ' tfa_filename_all '\n']);
        continue;
    else
        fprintf(['Start tfa for: ' tfa_filename_all '\n']);
    end

    %% -------------------------
    % Load data
    % -------------------------

    % Load data
    load(fullfile(prepro_path, prepro_filename));

    if exp_cfg.cutResponseWindow
        % Load log file
        behav = readtable(fullfile(beh_path, beh_filename), 'FileType', 'text');

        % Get sampling info
        fileNameIn = fullfile(exp_cfg.sourcedataPath, ['sub-', num2str(s)], ...
            'eeg', ['PEP_WP4_EEG', num2str(s), '.eeg']);
        hdr = ft_read_header(fileNameIn);
        Fs  = hdr.Fs;
        events = ft_read_event(fileNameIn);
        events = events(strcmp({events.type}, 'Stimulus'));
    end

    %% -------------------------
    % Loop blocks
    % -------------------------

    % Sore tfa of all blocks in one cell
    tf_all = cell(1, height(timelock.trial));

    for blk = 1:height(timelock.trial)

        % 5 bands: delta, theta, alpha, beta, gamma
        tf = cell(1, 5);

        % Adjust file name
        tfa_filename = ['PEP_WP4_EEG', num2str(s), '_tfa_block_',...
            num2str(blk), '.mat'];

        % Filter for current block
        cfg = [];
        cfg.trials = blk;
        data_block = ft_selectdata(cfg, timelock);

        % Remove NaN timepoints
        valid_time = data_block.time(~isnan(mean(data_block.trial, 2))); % check first channel
        cfg_trim = [];
        cfg_trim.latency = valid_time([1, end]);
        data_block = ft_selectdata(cfg_trim, data_block);

        % Get time points for tfa
        toi = -0.5:0.05:max(valid_time);

        %% =========================
        % DELTA (1–3 Hz)
        % =========================

        ft_cfg              = [];
        ft_cfg.trials       = 'all';
        ft_cfg.output       = 'pow';
        ft_cfg.pad          = 'nextpow2';
        ft_cfg.method       = 'mtmconvol';
        ft_cfg.keeptrials   = 'yes';
        ft_cfg.taper        = 'hanning';
        ft_cfg.foi          = 1:1:3;
        ft_cfg.toi          = toi;
        ft_cfg.t_ftimwin    = 5 ./ ft_cfg.foi;

        newTf = ft_freqanalysis(ft_cfg, data_block);

        if exp_cfg.normalize
            cfg = [];
            cfg.baseline     = [-0.5 -0.2];
            cfg.baselinetype = 'relchange';
            tf{1} = ft_freqbaseline(cfg, newTf);
        else
            tf{1} = newTf;
        end

        %% =========================
        % THETA (4–7 Hz)
        % =========================

        ft_cfg              = [];
        ft_cfg.trials       = 'all';
        ft_cfg.output       = 'pow';
        ft_cfg.pad          = 'nextpow2';
        ft_cfg.method       = 'mtmconvol';
        ft_cfg.keeptrials   = 'yes';
        ft_cfg.taper        = 'hanning';
        ft_cfg.foi          = 4:1:7;
        ft_cfg.toi          = toi;
        ft_cfg.t_ftimwin    = 5 ./ ft_cfg.foi;

        newTf = ft_freqanalysis(ft_cfg, data_block);

        if exp_cfg.normalize
            cfg = [];
            cfg.baseline     = [-0.5 -0.2];
            cfg.baselinetype = 'relchange';
            tf{2} = ft_freqbaseline(cfg, newTf);
        else
            tf{2} = newTf;
        end

        %% =========================
        % ALPHA (8–12 Hz)
        % =========================

        ft_cfg              = [];
        ft_cfg.trials       = 'all';
        ft_cfg.output       = 'pow';
        ft_cfg.pad          = 'nextpow2';
        ft_cfg.method       = 'mtmconvol';
        ft_cfg.keeptrials   = 'yes';
        ft_cfg.taper        = 'hanning';
        ft_cfg.foi          = 8:12;
        ft_cfg.toi          = toi;
        ft_cfg.t_ftimwin    = 5 ./ ft_cfg.foi;

        newTf = ft_freqanalysis(ft_cfg, data_block);

        if exp_cfg.normalize
            cfg = [];
            cfg.baseline     = [-0.5 -0.2];
            cfg.baselinetype = 'relchange';
            tf{3} = ft_freqbaseline(cfg, newTf);
        else
            tf{3} = newTf;
        end

        %% =========================
        % BETA (13–30 Hz)
        % =========================

        ft_cfg              = [];
        ft_cfg.trials       = 'all';
        ft_cfg.output       = 'pow';
        ft_cfg.pad          = 'nextpow2';
        ft_cfg.method       = 'mtmconvol';
        ft_cfg.keeptrials   = 'yes';
        ft_cfg.taper        = 'dpss';
        ft_cfg.foi          = 13:2:30;
        ft_cfg.toi          = toi;
        ft_cfg.t_ftimwin    = 5 ./ ft_cfg.foi;
        ft_cfg.tapsmofrq    = 0.4 .* ft_cfg.foi;

        newTf = ft_freqanalysis(ft_cfg, data_block);

        if exp_cfg.normalize
            cfg = [];
            cfg.baseline     = [-0.5 -0.2];
            cfg.baselinetype = 'relchange';
            tf{4} = ft_freqbaseline(cfg, newTf);
        else
            tf{4} = newTf;
        end

        %% =========================
        % GAMMA (31–70 Hz)
        % =========================

        ft_cfg              = [];
        ft_cfg.trials       = 'all';
        ft_cfg.output       = 'pow';
        ft_cfg.pad          = 'nextpow2';
        ft_cfg.method       = 'mtmconvol';
        ft_cfg.keeptrials   = 'yes';
        ft_cfg.taper        = 'dpss';
        ft_cfg.foi          = 31:2:70;
        ft_cfg.toi          = toi;
        ft_cfg.t_ftimwin    = 5 ./ ft_cfg.foi;
        ft_cfg.tapsmofrq    = 0.4 .* ft_cfg.foi;

        newTf = ft_freqanalysis(ft_cfg, data_block);

        if exp_cfg.normalize
            cfg = [];
            cfg.baseline     = [-0.5 -0.2];
            cfg.baselinetype = 'relchange';
            tf{5} = ft_freqbaseline(cfg, newTf);
        else
            tf{5} = newTf;
        end

        clear newTf

        %% -------------------------
        % Concatenate frequencies
        % -------------------------

        ft_cfg = [];
        ft_cfg.appenddim = 'freq';
        tf = ft_appendfreq(ft_cfg, tf{:});

        %% -------------------------
        % Save blocks tfa
        % -------------------------

        save(fullfile(tfa_path, tfa_filename), 'tf', '-v7.3');

        tf_all{blk} = tf;


        clear tf
    end
    %% -------------------------
    % Cut Response windows
    % -------------------------

    if exp_cfg.cutResponseWindow

        for blk = 1:height(timelock.trial)

            % Get current tf
            blkTf = tf_all{blk};

            % Get events of current block
            blkTrials = behav(behav.block == blk, :);
            blkEvents = events(blkTrials.trial + (blk - 1)*height(blkTrials));

            % Mark time windows that are good
            taskFreeTPs = true(1, length(blkTf.time));

            targetCount = 0;
            for i = 1:height(blkTrials)-1

                % Only trials WITH response
                if ~isnan(blkTrials.responseTime(i))
                    targetCount = targetCount +1;

                    % Define start/end in seconds
                    t_start = blkEvents(i).sample - blkEvents(1).sample + 250 + blkTrials.iti(i)*Fs; % + 500 ms
                    t_end   = blkEvents(i+1).sample - blkEvents(1).sample;

                    % Convert to samples
                    s_start = t_start / Fs;
                    s_end   = t_end   / Fs;

                    % Mark time points in this window
                    windowTPs = blkTf.time < s_end & blkTf.time > s_start;
                    taskFreeTPs(windowTPs) = false;
                end
            end
            disp(num2str(targetCount))

            % Filter for current block
            blkTf.time = blkTf.time(logical(taskFreeTPs));
            blkTf.trial = blk;
            blkTf.powspctrm = blkTf.powspctrm(:, :, :, logical(taskFreeTPs));
            tf_all{blk} = blkTf;

        end

        save(fullfile(tfa_path, tfa_filename_all), 'tf_all', '-v7.3');
    end

    clear timelock
end
end