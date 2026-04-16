function cfg = defaultCfg(cfg)

%% set default values for cfg

% exp parameters
if ~isfield(cfg, 'subNums'); cfg.subNums = 101:102;end
cfg.n = length(cfg.subNums);
if ~isfield(cfg, 'categories'); cfg.categories = {'bathroom', 'kitchen'};end
if ~isfield(cfg, 'exp_num'); cfg.exp_num = 1;end
if ~isfield(cfg, 'nTrials'); cfg.nTrials = 100;end

% eeg parameter
if ~isfield(cfg, 'nBlocks'); cfg.nBlocks = 20; end
if ~isfield(cfg, 'trialLength'); cfg.trialLength = 0.25;end
if ~isfield(cfg, 'ISC_types'); cfg.ISC_types = {'pairRep'};end
if ~isfield(cfg, 'channels'); cfg.channels = {'Fp1', 'F3', 'F7', 'FT9', 'FC5', 'FC1',...
        'C3', 'T7', 'TP9', 'CP5', 'CP1', 'Pz', 'P3', 'P7', 'O1', 'Oz', 'O2',...
        'P4', 'P8', 'TP10', 'CP6', 'CP2', 'Cz', 'C4', 'T8',...
        'FT10', 'FC6', 'FC2', 'F4', 'F8', 'Fp2', 'AF7', 'AF3', 'AFz',...
        'F1', 'F5', 'FT7', 'FC3', 'C1', 'C5', 'TP7', 'CP3', 'P1', 'P5', ...
        'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'P6', 'P2', 'CPz', 'CP4', 'TP8',...
        'C6', 'C2', 'FC4', 'FT8', 'F6', 'AF8', 'AF4', 'F2', 'FCz'};

% analsyis paramters
if ~isfield(cfg, 'dnn'); cfg.dnn = 'vgg16_imagenet';end
if ~isfield(cfg, 'analysis_names'); cfg.analysis_names = {'typical', 'control', 'photo'};end
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'saving'); cfg.saving = true; end
if ~isfield(cfg, 'dissimilarity'); cfg.dissimilarity = true; end
if ~isfield(cfg, 'classifier'); cfg.classifier = @cosmo_classify_lda; end

% define paths
if ~isfield(cfg, 'functionPath'); cfg.functionPath = fullfile(pwd,'utilities');end
if ~isfield(cfg, 'sourcedataPath'); cfg.sourcedataPath = fullfile(pwd, '..','sourcedata');end
if ~isfield(cfg, 'outputPath'); cfg.outputPath = fullfile(pwd, '..','derivatives');end
if ~isfield(cfg, 'figPath'); cfg.figPath = fullfile(pwd, '..', 'plots');end
if ~isfield(cfg, 'behPath'); cfg.behPath = fullfile(cfg.sourcedataPath, 'beh');end
if ~isfield(cfg, 'stimuliPath'); cfg.stimuliPath = fullfile(pwd, '..', 'stimuli');end
% cfg.sourcedataPath = fullfile(pwd, '..','sourcedata');
% cfg.outputPath = fullfile(pwd, '..','derivatives');
% cfg.figPath = fullfile(pwd, '..', 'plots');
% cfg.behPath = fullfile(cfg.sourcedataPath, 'beh');

% other standard configurations
if ~isfield(cfg, 'FontName'); cfg.FontName = 'Helvetica'; end
if ~isfield(cfg, 'FontSize'); cfg.FontSize = 15; end

% load colormap
load('utilities/colormaps.mat');
cfg.colormaps = colormaps;

end
