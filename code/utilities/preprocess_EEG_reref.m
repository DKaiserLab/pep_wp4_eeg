%% Housekeeping

close all
clear all
clc

%% Define Parameters
ss=[102 103 104 106 109 110 111 112 113 115 116 117 120 122];
prestim=0.5;
baseline=0.5;
poststim=0.5;

output_dir = fullfile(pwd, '..', '..', 'derivatives');

% enable fieldtrip functions
%restoredefaultpath;
%addpath ../../../../../../MATLAB/fieldtrip;
ft_defaults;

%% Start for-loop
for s=ss% for each subject    
           
    %% Specify input file

    fileName = fullfile(pwd, '..', '..', 'sourcedata', ['sub-', num2str(s)], 'eeg', ['PEP_WP4_EEG', num2str(s), '.eeg']);

    if exist(fullfile('..', '..', 'derivatives', ['sub-', num2str(s)], 'eeg', ['PEP_WP4_EEG', num2str(s), '_timelock_reref_w.mat']))==2
        continue
    end 

    %% Define Events

    cfg=[];
    cfg.dataset=fileName;
    cfg.trialdef.eventtype='Stimulus';
   
    cfg.trialdef.prestim=prestim;
    cfg.trialdef.poststim=poststim;
    cfg=ft_definetrial(cfg);
  
    %% Load Data and Preprocess Data
    
    cfg.hpfilter='no';
    cfg.lpfilter='no';
    cfg.bsfilter='yes';
    cfg.bsfreq=[48 52];

    % filter
%     cfg.bpfilter = 'yes';
%     cfg.bpfreq = [0.5 100];
%     cfg.bpfiltord = 4;
%     cfg.bpfiltype = 'but';

    cfg.demean='yes';
    cfg.baselinewindow=[-baseline,0];
    cfg.channel={'all', '-photo'};
    data=ft_preprocessing(cfg);

    %% Resample data 

    cfg = [];
    cfg.reref       = 'yes';
    cfg.reference_method = 'avg';
    cfg.refchannel  = 'all';          % Use all channels to compute average
    cfg.channel     = 'EEG';         
    data = ft_preprocessing(cfg, data);

    cfg=[];
    cfg.resamplefs=200;
    data=ft_resampledata(cfg,data);
    
    %% Look at the data to identify problematic trials
    
    cfg=[];
    cfg.viewmode='vertical';
    cfg=ft_databrowser(cfg,data);

    if isfield(cfg, 'trl')
        cfg = rmfield(cfg, 'trl');
    end

    data=ft_rejectartifact(cfg,data);

    %% Remove bad trials / channels
    
    cfg=[];
    cfg.showlabel='yes';
    cfg.method='summary'; 
    cfg.keepchannel='no';
    data=ft_rejectvisual(cfg,data);
    
    %% ICA     
    cfg = [];
    cfg.method='runica'; 
    cfg.numcomponent=30; % oder mehr?
    comp = ft_componentanalysis(cfg, data);

    % component topoplots
    figure;
    cfg=[];
    cfg.component = 1:30;%'all';
    layout = 'easycap-M1.txt';
    cfg.layout=layout; 
    cfg.comment='auto';%'no'
    ft_topoplotIC(cfg, comp);

    % component timecourse plots
    figure
    cfg=[];
    cfg.layout=layout; 
    cfg.viewmode='component';
    ft_databrowser(cfg,comp);

    % remove the bad components and backproject the data
    cfg = [];
    cfg.component = input('Which components do you want to remove? '); 
    cfg.demean='no';
    data=ft_rejectcomponent(cfg, comp, data);
    
    % define output_path, create folders
    output_path = fullfile(output_dir, ['sub-', num2str(s)], 'eeg');
    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end

    % transform to "timelocked" data and save the output
    cfg=[];
    cfg.outputfile = fullfile(output_path, ['PEP_WP4_EEG', num2str(s), '_timelock_reref_w']);
    cfg.keeptrials='yes';
    save(fullfile(output_dir, ['sub-', num2str(s)], 'eeg', ['PEP_WP4_EEG', num2str(s), '_timelock_reref_w']),'data');
    data=ft_timelockanalysis(cfg,data);
    
    % if you want: display the data
    %ft_databrowser(cfg,data);
    
end% subjects