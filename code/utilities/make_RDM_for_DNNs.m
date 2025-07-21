function d = make_RDM_for_DNNs(d, cfg, category, all_net_act)

% preparation
if ~isfield(cfg, 'plot_rdm'); cfg.plot_rdm = true; end
if ~isfield(cfg, 'dissimilarity'); cfg.dissimilarity = true; end

% get lyer name with _ instead of -
layer_name = strrep(cfg.layer_name, '-', '_');

if nargin < 4 % check how many input var
    % get feature map
    %     all_net_act = {d.DNN.(cfg.dnn).(cfg.analysis_name).(category).all_net_act.(layer_name).net_act};

    if ismember(cfg.dnn, {'dino','clip'})
        features_load_name = [char(cfg.analysis_name), '_', char(category), '.mat'];
    else
        features_load_name = [char(cfg.analysis_name), '_', char(category), '_', char(layer_name), '.mat'];
    end
    features_load_path = fullfile(pwd, '..', 'dnn_features', cfg.dnn, ['exp_', num2str(cfg.exp_num)], features_load_name);
    load(features_load_path)
    all_net_act = {dnn_features.all_net_act.net_act};
    all_net_act = cell2mat(all_net_act);

    % get images names
    all_img_names = {dnn_features.all_net_act.image_name};

end

% check if table
if istable(all_net_act)
    % get images names
    all_img_names = all_net_act.Properties.VariableNames;
end


%compute rdm
if cfg.plot_rdm == 1
    cfg.labels = [];
    cfg.plotting = 1;
    [fig, mat_out, ~] = make_RDM(all_net_act, cfg);
    fig_name = ['Correlation_', category, '_', cfg.analysis_name, '_layer_', layer_name, '_all_images'];
    title(strrep(fig_name, '_', ' '))

    % saveing figure
    fig_path = fullfile(cfg.figPath, ['exp_', num2str(cfg.exp_num)], cfg.analysis_name, category);
    if ~exist(fig_path, 'dir')
        mkdir(fig_path);
    end
    save_plot(fig, fig_name, fig_path)
else
    cfg.labels = [];
    cfg.plotting = 0;
    [~, mat_out, ~] = make_RDM(all_net_act, cfg);
end

% store RDMs for all image repeats in d
d.DNN.(cfg.dnn).(cfg.analysis_name).(category).all_images(cfg.layer).name = cfg.layer_name;
d.DNN.(cfg.dnn).(cfg.analysis_name).(category).all_images(cfg.layer).color = [0, 0, 0];
d.DNN.(cfg.dnn).(cfg.analysis_name).(category).all_images(cfg.layer).RDM = mat_out;


%% average for on subject level
mean_r = zeros(cfg.n,cfg.n);
for row = 1:cfg.n
    sub1 = cfg.subNums(row);

    for col = 1:cfg.n
        sub2 = cfg.subNums(col);

        sub1_indices = contains(all_img_names, sprintf('%0.3d', sub1));
        sub2_indices = contains(all_img_names, sprintf('%0.3d', sub2));
        cell_indices = (sub1_indices' & sub2_indices);
        mean_r(row,col) = mean(mean(mat_out(cell_indices)));

    end
end


if cfg.plot_rdm == 1
    % load colormap
    colormaps = load('white_zero_colormap.mat');

    % make RDM
    if cfg.dissimilarity
        if ~isfield(cfg, 'MinColorValue'); cfg.MinColorValue = 0; end
        if ~isfield(cfg, 'MaxColorValue'); cfg.MaxColorValue = 2; end
    else
        if ~isfield(cfg, 'MinColorValue'); cfg.MinColorValue = -1; end
        if ~isfield(cfg, 'MaxColorValue'); cfg.MaxColorValue = 1; end
    end

    % plot
    fig = figure;
    deoras_heatmap(mean_r, [], [], [],...
        'Colormap', colormaps.white_zero, 'Colorbar', true,...
        'MinColorValue', cfg.MinColorValue, 'MaxColorValue', cfg.MaxColorValue,...
        'TickAngle', 45, 'ShowAllTicks', 0);

    % give title
    fig_name = ['Correlation_', category, '_', cfg.analysis_name, '_layer_', cfg.layer_name, '_', cfg.dnn];
    title(strrep(fig_name, '_', ' '))

    % saveing figure
    fig_path = fullfile(cfg.figPath, ['exp_', num2str(cfg.exp_num)], cfg.analysis_name, category);
    if ~exist(fig_path, 'dir')
        mkdir(fig_path);
    end
    save_plot(fig, fig_name, fig_path)
end

% store RDMs for all image repeats in d
d.DNN.(cfg.dnn).(cfg.analysis_name).(category).subject_mean(cfg.layer).name = cfg.layer_name;
d.DNN.(cfg.dnn).(cfg.analysis_name).(category).subject_mean(cfg.layer).color = [0, 0, 0];
d.DNN.(cfg.dnn).(cfg.analysis_name).(category).subject_mean(cfg.layer).RDM = mean_r;
d.DNN.(cfg.dnn).(cfg.analysis_name).(category).subjects(cfg.layer) =  {cfg.subNums};


% %% average for on subject level
% mean_r = zeros(length(cfg.subNums_included),length(cfg.subNums_included));
% for row = 1:length(cfg.subNums_included)
%     sub1 = cfg.subNums_included(row);
% 
%     for col = 1:length(cfg.subNums_included)
%         sub2 = cfg.subNums_included(col);
% 
%         sub1_indices = contains(all_img_names, sprintf('%0.3d', sub1));
%         sub2_indices = contains(all_img_names, sprintf('%0.3d', sub2));
%         cell_indices = (sub1_indices' & sub2_indices);
%         mean_r(row,col) = mean(mean(mat_out(cell_indices)));
% 
%     end
% end
% 
% 
% if cfg.plot_rdm == 1
%     % load colormap
%     colormaps = load('white_zero_colormap.mat');
% 
%     % make RDM
%     if cfg.dissimilarity
%         if ~isfield(cfg, 'MinColorValue'); cfg.MinColorValue = 0; end
%         if ~isfield(cfg, 'MaxColorValue'); cfg.MaxColorValue = 2; end
%     else
%         if ~isfield(cfg, 'MinColorValue'); cfg.MinColorValue = -1; end
%         if ~isfield(cfg, 'MaxColorValue'); cfg.MaxColorValue = 1; end
%     end
% 
%     % plot
%     fig=figure;
%     deoras_heatmap(mean_r, [], [], [],...
%         'Colormap', colormaps.white_zero, 'Colorbar', true,...
%         'MinColorValue', cfg.MinColorValue, 'MaxColorValue', cfg.MaxColorValue,...
%         'TickAngle', 45, 'ShowAllTicks', 0);
% 
%     % give title
%     fig_name = ['Correlation_', category, '_', cfg.analysis_name, '_layer_', cfg.layer_name, '_', cfg.dnn];
%     title(strrep(fig_name, '_', ' '))
% 
%     % saveing figure
%     fig_path = fullfile(cfg.figPath, ['exp_', num2str(cfg.exp_num)], cfg.analysis_name, category);
%     save_plot(fig, fig_name, fig_path)
% end
% 
% % store RDMs for all image repeats in d
% d.DNN.(cfg.dnn).(cfg.analysis_name).(category).subject_mean(cfg.layer).name = cfg.layer_name;
% d.DNN.(cfg.dnn).(cfg.analysis_name).(category).subject_mean(cfg.layer).color = [0, 0, 0];
% d.DNN.(cfg.dnn).(cfg.analysis_name).(category).subject_mean(cfg.layer).RDM = mean_r;
% d.DNN.(cfg.dnn).(cfg.analysis_name).(category).subjects(cfg.layer) =  {cfg.subNums};

end
