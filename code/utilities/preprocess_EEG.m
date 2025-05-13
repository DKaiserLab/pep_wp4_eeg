%% Housekeeping

close all
clear all
clc

%% Define Parameters

ss=109; %[106 109 111 113 110]
prestim=0.5;
baseline=0.5;
poststim=0.5;
<<<<<<< HEAD
=======
%output_dir='Users/theaschmitt/Documents/Data Science/6. Semester/Thesis/Test/';
%output_dir = fullfile('Users', 'theaschmitt', 'Documents', 'Data Science', '6. Semester', 'Thesis', 'Test');
>>>>>>> cf25af5959fdfd67a5ba0d34563c2610b6511b1b
output_dir = fullfile(pwd, '..', '..', 'derivatives');

% enable fieldtrip functions
%restoredefaultpath;
%addpath ../../../../../../MATLAB/fieldtrip;
ft_defaults;

%% Start for-loop
for s=ss% for each subject    
           
    %% Specify input file
    
<<<<<<< HEAD
    fileName = fullfile(pwd, '..', '..', 'sourcedata', ['sub-', num2str(s)], 'eeg', ['PEP_WP4_EEG', num2str(s), '.eeg']);
    
=======
    %fileName=['Users/theaschmitt/Documents/Data Science/6. Semester/Thesis/Daten/','PEP_WP4_EEG', num2str(s,'%02.f'),'.eeg'];
    fileName = fullfile(pwd, '..', '..', 'sourcedata', ['sub-', num2str(s)], 'eeg', ['PEP_WP4_EEG', num2str(s), '.eeg']);
>>>>>>> cf25af5959fdfd67a5ba0d34563c2610b6511b1b
    %% Define Events

    cfg=[];
    cfg.dataset=fileName;
    cfg.trialdef.eventtype='Stimulus';
   
    cfg.trialdef.prestim=prestim;
    cfg.trialdef.poststim=poststim;
    cfg=ft_definetrial(cfg);
%     
%     two_dig_logical=cfg.trl(:,4)<100;
%     two_dig_triggers=cfg.trl(two_dig_logical,4);
%     cfg.trl(two_dig_logical,4)=cellfun(@(c)c(1),""+two_dig_triggers)-'0';
%     three_dig_triggers=cfg.trl(~two_dig_logical,4);
%     cfg.trl(~two_dig_logical,4)=cellfun(@(c)c(2),""+three_dig_triggers)-'0';
%     
    %% Load Data and Preprocess Data
    
    cfg.hpfilter='no';
    cfg.lpfilter='no';
    cfg.bsfilter='yes';
    cfg.bsfreq=[48 52];
    cfg.demean='yes';
    cfg.baselinewindow=[-baseline,0];
    cfg.channel={'all', '-photo'};
    data=ft_preprocessing(cfg);

    %% Resample data 
    
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
    cfg.method='fastica'; 
    cfg.fastica.numOfIC=30;
    comp = ft_componentanalysis(cfg, data);

    % component topoplots
    figure
    cfg=[];
    cfg.component=1:30;   
    layout = 'easycap-M1.txt';
    cfg.layout=layout; 
    cfg.comment='no';
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
<<<<<<< HEAD
    cfg.outputfile = fullfile(output_path, ['PEP_WP4_EEG', num2str(s), '_timelock']);
=======
    cfg.outputfile=fullfile(output_dir, ['sub-', num2str(s)], 'eeg', ['PEP_WP4_EEG', num2str(s), '_timelock']);
>>>>>>> cf25af5959fdfd67a5ba0d34563c2610b6511b1b
    cfg.keeptrials='yes';
    save(fullfile(output_dir, ['sub-', num2str(s)], 'eeg', ['PEP_WP4_EEG', num2str(s), '_timelock']),'data');
    data=ft_timelockanalysis(cfg,data);
    
    % if you want: display the data
    %ft_databrowser(cfg,data);
    
end% subjects