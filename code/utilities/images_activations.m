function [d, dnn_features] = images_activations(cfg, d, images, dnn_features, net, layer_name, category)

ix=0;
for image_name = images.(cfg.analysis_name).(category).image_names
    %image counter
    ix=ix+1;

    % check if it already exists and skip is thats the case
    if exist('dnn_features', 'var') % if exist(output_dir, 'var')
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
    new_net_act = activations(net, img, cfg.layer_name); % compute activation of the current layer
    dnn_features.all_net_act(ix).image_name = char(image_name); % save image-name
    dnn_features.all_net_act(ix).net_act = new_net_act(:); % and activation for the current layer

    d.DNN.(cfg.dnn).(cfg.analysis_name).(category).all_net_act.(layer_name)(ix).image_name = char(image_name);
    d.DNN.(cfg.dnn).(cfg.analysis_name).(category).all_net_act.(layer_name)(ix).net_act = new_net_act(:);
end
end