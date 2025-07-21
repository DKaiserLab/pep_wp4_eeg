function [d, cfg] = get_dnn_layer_similarity(d, images, cfg)

% evaluate input
if ~isfield(cfg, 'RDM_of_DNNS'); cfg.RDM_of_DNNS = 1; end
if ~isfield(cfg, 'layer_type'); cfg.layer_type = 'late';end
if ~isfield(cfg, 'dnn'); cfg.dnn = 'vgg16_imagenet';end
if ~isfield(cfg, 'categories');  cfg.categories = {'kitchen', 'bathroom'};end
if ~isfield(cfg, 'analysis_names');  cfg.analysis_names = {'typical', 'control', 'photo'};end

% load network
[cfg, net] = load_dnn(cfg);

%% run through dnn layers
for layer=1:length(cfg.loi)
    cfg.layer = layer;

    %find current layer name
    cfg.layer_name=net.Layers(cfg.loi(layer)).Name;
    cfg.(cfg.dnn).all_layer_names = net.Layers.Name; % ?

    %cycle through categories
    for category = cfg.categories
        category = char(category);

        % loop through analysis
        for analysis_name = cfg.analysis_names
            cfg.analysis_name = char(analysis_name);

            % check if features exit
            layer_name = strrep(cfg.layer_name, '-', '_'); % customize name
            features_save_name = [char(cfg.analysis_name), '_', char(category), '_', char(layer_name), '.mat']; % mat-name
            features_save_path = fullfile(pwd, '..', 'dnn_features', cfg.dnn, ['exp_', num2str(cfg.exp_num)]);
            file_name = fullfile(features_save_path, features_save_name); % compose path 

            if ~exist(features_save_path, 'dir')
                mkdir(features_save_path)
            end

            % if feature activations exist load them
            if exist(file_name, 'file')
                load(file_name)

                % replace by filtered struct with only subjects that are included
                [cfg, dnn_features] = filter_existing_images(cfg, dnn_features);

            else
                dnn_features = struct;
                %                     for I=1:numel(images.(cfg.analysis_name).(category).image_names)
                %                         dnn_features.all_net_act(I).image_name = [];
                %                         dnn_features.all_net_act(I).net_act = [];
                %                     end

                % calculate and store activations of DNN for each image
                [d, dnn_features] = images_activations(cfg, d, images, dnn_features, net, layer_name, category);
            end
          
            %show progress
            disp(['analysis ', cfg.analysis_name,' - layer #',num2str(layer), ' - category: ', category])

            % save features
            save(file_name, 'dnn_features')

            if cfg.RDM_of_DNNS == 1
                % get RDMs for all images and for subject averaged images
                %[cfg, dnn_features] = filter_images(cfg, dnn_features);
                d = make_RDM_for_DNNs(d, cfg, category);  
            end
            
        end % analysis, added
    end % category
end % layer
end