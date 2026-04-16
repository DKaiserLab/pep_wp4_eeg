function preprocess_EEG(cfg)

%% Define Parameters
prestim=1;
baseline=0.5;
poststim=1;


% enable fieldtrip functions
%restoredefaultpath;
%addpath ../../../../../../MATLAB/fieldtrip;
ft_defaults;

%% Start for-loop
for s=cfg.subNums% for each subject

    %% Specify input and output file name

    fileNameIn = fullfile(cfg.sourcedataPath, ['sub-', num2str(s)], 'eeg', ['PEP_WP4_EEG', num2str(s), '.eeg']);
    fileNameOut = fullfile(cfg.outputPath, ['sub-', num2str(s)], 'eeg', ['PEP_WP4_EEG', num2str(s), '_timelock_reref_s2_filtered.mat']);

    if exist(fileNameOut, 'file')
        disp(['Pre-processed data for ', ['sub-', num2str(s)], ' exist already', newline]);
        continue
    end

    %% Define Events

    cfg_temp=[];
    cfg_temp.dataset=fileNameIn;
    cfg_temp.trialdef.eventtype='Stimulus';

    cfg_temp.trialdef.prestim=prestim;
    cfg_temp.trialdef.poststim=poststim;
    cfg_temp=ft_definetrial(cfg_temp);


    disp(['Start preprocessing for sub-', num2str(s)]);
    
    %% Load Data and Preprocess Data

    %cfg_temp.hpfilter='no';
    %cfg_temp.lpfilter='no';
    cfg_temp.bsfilter='yes';
    cfg_temp.bsfreq=[48 52];

    % filtering
    cfg_temp.hpfilter = 'yes';
    cfg_temp.hpfreq   = 0.5;
    cfg_temp.hpfilttype  = 'firws';

    cfg_temp.lpfilter = 'yes';
    cfg_temp.lpfreq   = 100;
    cfg_temp.hpfilttype  = 'firws';

    % filter
    %     cfg_temp.bpfilter = 'yes';
    %     cfg_temp.bpfreq = [0.5 100];
    %     cfg_temp.bpfiltord = 4;
    %     cfg_temp.bpfiltype = 'but';

    cfg_temp.demean='yes';
    cfg_temp.baselinewindow=[-baseline,0];
    cfg_temp.channel={'all', '-photo'};
    data=ft_preprocessing(cfg_temp);

    %% Reref
    cfg_temp = [];
    cfg_temp.reref       = 'yes';
    cfg_temp.reference_method = 'avg';
    cfg_temp.refchannel  = 'all';          % Use all channels to compute average
    cfg_temp.channel     = 'EEG';
    data = ft_preprocessing(cfg_temp, data);
    
    %% Resample data
    cfg_temp=[];
    cfg_temp.resamplefs=200;
    data=ft_resampledata(cfg_temp,data);

    %% Look at the data to identify problematic trials

    cfg_temp=[];
    cfg_temp.viewmode='vertical';
    cfg_temp=ft_databrowser(cfg_temp,data);

    if isfield(cfg_temp, 'trl')
        cfg_temp = rmfield(cfg_temp, 'trl');
    end

    data=ft_rejectartifact(cfg_temp,data);

    %% Remove bad trials / channels

    cfg_temp=[];
    cfg_temp.showlabel='yes';
    cfg_temp.method='summary';
    cfg_temp.keepchannel='no';
    data=ft_rejectvisual(cfg_temp,data);

    %% ICA
    cfg_temp = [];
    cfg_temp.method='runica';
    cfg_temp.numcomponent=30;
    comp = ft_componentanalysis(cfg_temp, data);

    % component topoplots
    figure;
    cfg_temp=[];
    cfg_temp.component = 1:30; %'all';
    layout = 'easycap-M1.txt';
    cfg_temp.layout=layout;
    cfg_temp.comment='auto'; %'no'
    ft_topoplotIC(cfg_temp, comp);

    % component timecourse plots
    figure
    cfg_temp=[];
    cfg_temp.layout=layout;
    cfg_temp.viewmode='component';
    ft_databrowser(cfg_temp,comp);

    % remove the bad components and backproject the data
    cfg_temp = [];
    cfg_temp.component = input('Which components do you want to remove? ');
    cfg_temp.demean='no';
    data=ft_rejectcomponent(cfg_temp, comp, data);

    %% define output_path, create folders
    output_path = fullfile(cfg.outputPath, ['sub-', num2str(s)], 'eeg');
    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end

    %% transform to "timelocked" data and save the output
    cfg_temp=[];
    cfg_temp.outputfile = fileNameOut;
    cfg_temp.keeptrials='yes';
    timelock = ft_timelockanalysis(cfg_temp, data);
    save(fileNameOut,'timelock', "comp");
    
    % if you want: display the data
    %ft_databrowser(cfg,data);

end% subjects
end