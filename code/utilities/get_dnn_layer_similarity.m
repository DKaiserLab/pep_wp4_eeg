function [d, cfg] = get_dnn_layer_similarity(d, images, cfg)

% evaluate input
if ~isfield(cfg, 'RDM_of_DNNS'); cfg.RDM_of_DNNS = 1; end
if ~isfield(cfg, 'layer_type'); cfg.layer_type = 'late';end
if ~isfield(cfg, 'dnn'); cfg.dnn = 'vgg16_imagenet';end
if ~isfield(cfg, 'categories');  cfg.categories = {'kitchen', 'bathroom'};end
if ~isfield(cfg, 'analysis_names');  cfg.analysis_names = {'typical', 'control'};end
dnn_features = struct;

% load network
if strcmp(cfg.dnn,'googlenet_places365')
    net=googlenet('weights','places365');

    % get layers
    if strcmp(cfg.layer_type,'early_mid_late')
        cfg.loi = [8,68,139];
    elseif strcmp(cfg.layer_type,'late')
        cfg.loi= 139;
    elseif strcmp(cfg.layer_type,'mid_late')
        cfg.loi= [68,139];
    elseif strcmp(cfg.layer_type,'all_output')
        cfg.loi=[8,25,39,54,68,82,96,110,125,139,142];
    elseif strcmp(cfg.layer_type,'all')
        cfg.loi = 1:144;
    end

elseif contains(cfg.dnn,'vgg')

    if strcmp(cfg.dnn,'vggnet16_places365')
        load(fullfile('..', 'vgg16_places365/vggnet16_places365.mat'));

    elseif strcmp(cfg.dnn,'vgg16_imagenet')
        load(fullfile('..', 'vgg16_imagenet/vgg16.mat'));
        net=vgg16;
        clear vgg16
    end

    % get layers
    if strcmp(cfg.layer_type,'early_mid_late')
        cfg.loi = [4,21,36];
    elseif strcmp(cfg.layer_type,'late')
        cfg.loi = 36;
    elseif strcmp(cfg.layer_type,'mid_late')
        cfg.loi = [21,36];
    elseif strcmp(cfg.layer_type,'all_output')
        %Get the Conv and FC Layers
        lx=0;
        for layer=1:length(net.Layers)
            if strfind(net.Layers(layer).Name,'conv')==1
                lx=lx+1;
                cfg.loi(lx)=layer;
            elseif strfind(net.Layers(layer).Name,'fc')==1
                lx=lx+1;
                cfg.loi(lx)=layer;
            end
        end
    elseif strcmp(cfg.layer_type,'all')
        cfg.loi = 1:41;
    end
end

%run through dnn layers
for layer=1:length(cfg.loi)
    cfg.layer = layer;

    %find current layer name
    cfg.layer_name=net.Layers(cfg.loi(layer)).Name;
    cfg.(cfg.dnn).all_layer_names =net.Layers.Name;

    %cycle through categories
    for category = cfg.categories
        category = char(category);

        % check if features exit 
        layer_name = strrep(cfg.layer_name, '-', '_');
        features_save_name = [char(cfg.analysis_name), '_', char(category), '_', char(layer_name), '.mat'];
        features_save_path = fullfile(pwd, '..', 'dnn_features', cfg.dnn, ['exp_', num2str(cfg.exp_num)]);
        file_name = fullfile(features_save_path, features_save_name);
        if ~exist(features_save_path, 'dir')
            mkdir(features_save_path)
        else
            if exist(file_name, 'file')
                load(file_name)

                % Loop through struct and extract machting subject numbers
                count = 0;
                for i = 1:length(dnn_features.all_net_act)

                    % Get the current image name and subject number
                    current_img_name = dnn_features.all_net_act(i).image_name;
                    number = regexp(current_img_name, '\d+', 'match');
                    number = str2double(number{1});

                    % Check if the number is part of the specified array
                    if ismember(number, cfg.subNums)
                        count = count + 1;
                        
                        % If yes, add the field to the filtered structure
                        filtered_struct.all_net_act(count).image_name = current_img_name;
                        filtered_struct.all_net_act(count).net_act =...
                            dnn_features.all_net_act(i).net_act;
                    end
                end

                % replace by filtered struct with only subject that are
                % included
                dnn_features = filtered_struct;
            else
                dnn_features = struct;
                for I=1:numel(images.(cfg.analysis_name).(category).image_names)
                    dnn_features.all_net_act(I).image_name = [];
                    dnn_features.all_net_act(I).net_act = [];
                end
            end 
        end 

        %show progress
        disp(['analysis ', cfg.analysis_name,' - layer #',num2str(layer), ' - category: ', category])

        %image counter
        ix=0;

        for image_name = images.(cfg.analysis_name).(category).image_names

            %image counter
            ix=ix+1;

            % check if it already exists and skip is thats the case
            if exist('dnn_features', 'var')
                if isfield(dnn_features, 'all_net_act')
                    if isfield(dnn_features.all_net_act, 'image_name')
                        if ix <= length({dnn_features.all_net_act.image_name})
                            if strcmp(image_name, dnn_features.all_net_act(ix).image_name)
                                continue
                            end
                        end
                    end
                end
            end

            % get image
            img = images.(cfg.analysis_name).(category).image_data{ix};

            %preprocess image
            img=imresize(img,[227,227]); % resize the image
            if size(img,3)==1
                img=cat(3,img,img,img); % adjust image dimensions
            end

            %calculate and store activations
            new_net_act = activations(net, img, cfg.layer_name);
            dnn_features.all_net_act(ix).image_name = char(image_name);
            dnn_features.all_net_act(ix).net_act = new_net_act(:);

            %             d.DNN.(cfg.dnn).(cfg.analysis_name).(category).all_net_act.(layer_name)(ix).image_name = char(image_name);
            %             d.DNN.(cfg.dnn).(cfg.analysis_name).(category).all_net_act.(layer_name)(ix).net_act = new_net_act(:);

        end

        % save features
        save(file_name, 'dnn_features')

        if cfg.RDM_of_DNNS == 1
            % get RDMs for all images and for subject averaged images
            d = make_RDM_for_DNNs(d, cfg, category);
        end
    end % category
end % layer
end