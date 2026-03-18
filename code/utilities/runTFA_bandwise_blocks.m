function runTFA_bandwise_blocks(exp_cfg)

if ~isfield(exp_cfg, 'normalize'); exp_cfg.normalize = false; end

for s = exp_cfg.subNums

    %% -------------------------
    % Paths
    % -------------------------

    prepro_path = fullfile(exp_cfg.outputPath, ['sub-', num2str(s)], 'eeg');
    tfa_path = prepro_path;
    prepro_filename = ['PEP_WP4_EEG', num2str(s), '_blocks.mat'];
    tfa_filename = ['PEP_WP4_EEG', num2str(s), '_tfa_block_',...
        num2str(1), '.mat'];

    %% -------------------------
    % Skip if exists
    % -------------------------

    if exist(fullfile(tfa_path, tfa_filename), 'file')
        fprintf(['File already exists: ' tfa_filename '\n']);
        continue;
    else
        fprintf(['Start tfa for: ' tfa_filename '\n']);
    end

    %% -------------------------
    % Load data
    % -------------------------

    load(fullfile(prepro_path, prepro_filename));

    for iBlk = 1:height(timelock.trial)

        % Now 5 bands: delta, theta, alpha, beta, gamma
        tf = cell(1, 5);

        % Adjust file name
        tfa_filename = ['PEP_WP4_EEG', num2str(s), '_tfa_block_',...
            num2str(iBlk), '.mat'];

        % Filter for current block
        cfg = [];
        cfg.trials = iBlk;
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

        %% -------------------------
        % Cleanup
        % -------------------------

        clear newTf

        %% -------------------------
        % Concatenate frequencies
        % -------------------------

        ft_cfg = [];
        ft_cfg.appenddim = 'freq';
        tf = ft_appendfreq(ft_cfg, tf{:});

        %% -------------------------
        % Save
        % -------------------------

        save(fullfile(tfa_path, tfa_filename), 'tf', '-v7.3');

        clear tf

    end

    clear timelock
end
end