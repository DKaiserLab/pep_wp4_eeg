function preprocess_EEG_blocks(cfg)

%% =========================
%  EEG Block-wise Preprocessing
%  =========================

%% Define Parameters

prestim  = 1;    
poststim = 2;    

ft_defaults;

%% =========================
%  Loop over subjects
% =========================

for s = cfg.subNums

    %% -------------------------
    % File paths
    % -------------------------
    
    fileNameIn = fullfile(cfg.sourcedataPath, ['sub-', num2str(s)], ...
        'eeg', ['PEP_WP4_EEG', num2str(s), '.eeg']);
    
    fileNameOut = fullfile(cfg.outputPath, ['sub-', num2str(s)], ...
        'eeg', ['PEP_WP4_EEG', num2str(s), '_blocks.mat']);

    if exist(fileNameOut, 'file')
        disp(['Block data for sub-', num2str(s), ' already exists', newline]);
        continue
    end

    disp(['Start BLOCK preprocessing for sub-', num2str(s)]);

    %% -------------------------
    % Load cleaned trial-wise data → get good channels
    % -------------------------
    
    clean_file = fullfile(cfg.outputPath, ['sub-', num2str(s)], ...
        'eeg', ['PEP_WP4_EEG', num2str(s), '_timelock_reref_s2.mat']);
    
    load(clean_file, 'timelock');
    good_channels = timelock.label;

    %% -------------------------
    % Define BLOCK trials
    % -------------------------
    
    hdr   = ft_read_header(fileNameIn);
    event = ft_read_event(fileNameIn);

    stim_events = event(strcmp({event.type}, 'Stimulus'));
    samples = [stim_events.sample];

    trials_per_block = 100;
    n_blocks = floor(length(samples) / trials_per_block);

    trl = [];

    for b = 1:n_blocks
        
        start_sample = samples((b-1)*trials_per_block + 1);
        end_sample   = samples(b*trials_per_block);
        
        prestim_samp  = round(prestim  * hdr.Fs);
        poststim_samp = round(poststim * hdr.Fs);
        
        trlbegin = start_sample - prestim_samp;
        trlend   = end_sample   + poststim_samp;
        offset   = -prestim_samp;
        
        trl = [trl; trlbegin trlend offset];
    end

    %% -------------------------
    % Preprocessing (MAIN DATA)
    % -------------------------
    
    cfg_temp = [];
    cfg_temp.dataset = fileNameIn;
    cfg_temp.trl     = trl;

    % filtering
    cfg_temp.hpfilter = 'yes';
    cfg_temp.hpfreq   = 0.5;
    cfg_temp.hpfilttype  = 'firws';

    cfg_temp.lpfilter = 'yes';
    cfg_temp.lpfreq   = 100;
    cfg_temp.hpfilttype  = 'firws';

    % Notch
    cfg_temp.bsfilter='yes';
    cfg_temp.bsfreq=[48 52];

    % Demean (no trial baseline)
    cfg_temp.demean = 'yes';
    cfg_temp.baselinewindow = 'all';

    % Use same cleaned channels
    cfg_temp.channel = good_channels;

    data = ft_preprocessing(cfg_temp);

    %% -------------------------
    % Re-reference
    % -------------------------
    
    cfg_temp = [];
    cfg_temp.reref       = 'yes';
    cfg_temp.reference_method = 'avg';
    cfg_temp.refchannel  = 'all';
    cfg_temp.channel     = 'EEG';

    data = ft_preprocessing(cfg_temp, data);

    %% -------------------------
    % Resample
    % -------------------------
    
    cfg_temp = [];
    cfg_temp.resamplefs = 200;

    data = ft_resampledata(cfg_temp, data);

    %% -------------------------
    % Artifact rejection (visual)
    % -------------------------
    
    cfg_temp = [];
    cfg_temp.viewmode = 'vertical';
    % cfg_temp = ft_databrowser(cfg_temp, data);

    if isfield(cfg_temp, 'trl')
        cfg_temp = rmfield(cfg_temp, 'trl');
    end

    data = ft_rejectartifact(cfg_temp, data);

    %% -------------------------
    % ICA (on filtered copy!)
    % -------------------------
    
    % Create ICA-optimized copy
    cfg_ica = [];
    cfg_ica.hpfilter = 'yes';
    cfg_ica.hpfreq   = 1;   % stronger HP for ICA

    data_ica = ft_preprocessing(cfg_ica, data);

    % Run ICA
    cfg_temp = [];
    cfg_temp.method = 'runica';
    cfg_temp.numcomponent = 30;

    comp = ft_componentanalysis(cfg_temp, data_ica);

    %% -------------------------
    % Inspect components
    % -------------------------
    
    figure;
    cfg_temp = [];
    cfg_temp.component = 1:30;
    cfg_temp.layout = 'easycap-M1.txt';
    cfg_temp.comment = 'auto';
    ft_topoplotIC(cfg_temp, comp);

    figure;
    cfg_temp = [];
    cfg_temp.layout = 'easycap-M1.txt';
    cfg_temp.viewmode = 'component';
    ft_databrowser(cfg_temp, comp);

    %% -------------------------
    % Remove components (applied to ORIGINAL data!)
    % -------------------------
    
    cfg_temp = [];
    cfg_temp.component = input('Which components do you want to remove? ');
    cfg_temp.demean = 'no';

    data = ft_rejectcomponent(cfg_temp, comp, data);

    %% -------------------------
    % Save data
    % -------------------------
    
    output_path = fullfile(cfg.outputPath, ['sub-', num2str(s)], 'eeg');
    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end

    cfg_temp = [];
    cfg_temp.outputfile = fileNameOut;
    cfg_temp.keeptrials = 'yes';

    timelock = ft_timelockanalysis(cfg_temp, data);

    save(fileNameOut, 'timelock');

end

end