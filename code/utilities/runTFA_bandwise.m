function runTFA_bandwise(exp_cfg)

if ~isfield(exp_cfg, 'normalize'); exp_cfg.normalize = false; end


% for s = nSubjs
for s = exp_cfg.subNums


    % get paths and filenames
    prepro_path = fullfile(exp_cfg.outputPath, ['sub-', num2str(s)], 'eeg');
    prepro_filename = ['PEP_WP4_EEG', num2str(s), '_timelock_reref_s2.mat'];
    tfa_path = fullfile(exp_cfg.outputPath, ['sub-', num2str(s)], 'eeg');
    tfa_filename = ['PEP_WP4_EEG', num2str(s), '_tfa.mat'];


    tf = cell(1, 3); % Initialize the cell array to store the results

    % Check if the file already exists - if it does, skip it
    if exist(fullfile(tfa_path, tfa_filename), 'file')
        fprintf(['File already exists: ' tfa_filename '\n']);
        continue;
    else
        fprintf(['Start tfa for: ' tfa_filename '\n']);
    end

    %% eeg timelock
    load(fullfile(prepro_path, prepro_filename));

    % select timelock
    % ft_cfg              = [];
    % ft_cfg.latency      = [0, 3]; % 0-3s
    % timelock = ft_selecttimelock(ft_cfg, timelock);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % alpha band
    ft_cfg              = [];
    ft_cfg.trials       = 'all';
    ft_cfg.output       = 'pow';
    ft_cfg.pad          = 'nextpow2';
    ft_cfg.method       = 'mtmconvol';
    ft_cfg.keeptrials   = 'yes';
    ft_cfg.taper        = 'hanning'; % single taper
    ft_cfg.foi          = 8:12;
    ft_cfg.toi          = -0.5:0.02:0.5;
    ft_cfg.t_ftimwin    = 4 ./ ft_cfg.foi;
    newTf = ft_freqanalysis(ft_cfg, timelock);

    if exp_cfg.normalize
        cfg = [];
        cfg.baseline     = [-0.5 -0.2];
        cfg.baselinetype = 'relchange';
        tf{1} = ft_freqbaseline(cfg, newTf);
    else
        tf{1} = newTf;
    end

    % beta band
    ft_cfg              = [];
    ft_cfg.trials       = 'all';
    ft_cfg.output       = 'pow';
    ft_cfg.pad          = 'nextpow2';
    ft_cfg.method       = 'mtmconvol';
    ft_cfg.keeptrials   = 'yes';
    ft_cfg.taper        = 'dpss';
    ft_cfg.foi          = 13:2:30;
    ft_cfg.toi          = -0.5:0.02:0.5;
    ft_cfg.t_ftimwin    = 5 ./ ft_cfg.foi;
    ft_cfg.tapsmofrq    = 0.4 .* ft_cfg.foi;
    newTf = ft_freqanalysis(ft_cfg, timelock);

    if exp_cfg.normalize
        cfg = [];
        cfg.baseline     = [-0.5 -0.2];
        cfg.baselinetype = 'relchange';
        tf{2} = ft_freqbaseline(cfg, newTf);
    else
        tf{2} = newTf;
    end

    % gamma band
    ft_cfg              = [];
    ft_cfg.trials       = 'all';
    ft_cfg.output       = 'pow';
    ft_cfg.pad          = 'nextpow2';
    ft_cfg.method       = 'mtmconvol';
    ft_cfg.keeptrials   = 'yes';
    ft_cfg.taper        = 'dpss'; % multiple tapers
    ft_cfg.foi          = 31:2:70;
    ft_cfg.tapsmofrq    = 0.4 .* ft_cfg.foi;
    ft_cfg.toi          = -0.5:0.02:0.5;
    ft_cfg.t_ftimwin    = 6 ./ ft_cfg.foi;
    ft_cfg.tapsmofrq    = 0.4 .* ft_cfg.foi;
    newTf = ft_freqanalysis(ft_cfg, timelock);

    if exp_cfg.normalize
        cfg = [];
        cfg.baseline     = [-0.5 -0.2];
        cfg.baselinetype = 'relchange';
        tf{3} = ft_freqbaseline(cfg, newTf);
    else
        tf{3} = newTf;
    end

    % free memory
    clear timelock_sel
    clear timelock
    clear newTf

    % concatenate data
    ft_cfg = [];
    ft_cfg.appenddim = 'freq';
    tf = ft_appendfreq(ft_cfg, tf{:});

    % save data
    save(fullfile(tfa_path, tfa_filename), 'tf', '-v7.3');

    clear tf

end
end