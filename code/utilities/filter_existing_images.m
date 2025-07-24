function [cfg, filtered_struct] = filter_existing_images(cfg, dnn_features)
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
        filtered_struct.all_net_act(count).net_act = dnn_features.all_net_act(i).net_act;
        cfg.subNums_included(i)=number;
    end
end

cfg.subNums_included = unique(cfg.subNums_included);
cfg.subNums_included(cfg.subNums_included == 0) = [];

% check if all subjects have feature activations
for iSub = cfg.subNums
    if ~ismember(iSub, cfg.subNums_included)
        warning(['Missing feature activations for ', num2str(iSub)])
    end
end
end