function [cfg, d, images] = make_RDM_for_stimuli(cfg, d, images)

%% load stimuli
stimuli_image_folder = fullfile(cfg.stimuliPath);

if ~isfield(images, 'stimuli')
    images.stimuli = struct();
end

for category = cfg.categories
    category = char(category);

    files = dir(fullfile(stimuli_image_folder, category, '*'));
    files = files(~[files.isdir]);
    files = struct2table(files);

    % save
    filelist.(category) = files;

    for i = 1:height(filelist.(category))

        % get image name
        field_name = char(filelist.(category).name(i));

        % load image
        image_path = fullfile(filelist.(category).folder(i), field_name);
        new_image = imread(char(image_path));

        % make gray scale
        if ndims(new_image) == 3
            if size(new_image,3) == 4
                new_image = rgb2gray(new_image(:,:,1:3));
            else
                new_image = rgb2gray(new_image);
            end
        end

        % Store the image data in structre
        images.stimuli.(category).image_data{i} = new_image;
        images.stimuli.(category).image_names{i} = field_name;

        % show progress
        disp(['Loaded image: ',field_name])
    end

    
end

%% get DNN layer activations
[d, cfg] = get_dnn_layer_similarity(d, images, cfg);
% d.DNN.vgg16_imagenet.stimuli.all = struct;
% for i=1:length(cfg.loi)
%     layer = cfg.layer_names{i};
%     d.DNN.vgg16_imagenet.stimuli.all.all_images.(layer).RDM = nan(100,100);
%     d.DNN.vgg16_imagenet.stimuli.all.all_images.(layer).RDM(1:50, 1:50) = d.DNN.vgg16_imagenet.stimuli.bathroom.all_images(i).RDM;
%     d.DNN.vgg16_imagenet.stimuli.all.all_images.(layer).RDM(51:100, 51:100) = d.DNN.vgg16_imagenet.stimuli.kitchen.all_images(i).RDM;
%     d.DNN.vgg16_imagenet.stimuli.all.all_images.(layer).name = ['stimuli_',layer];
%     d.DNN.vgg16_imagenet.stimuli.all.all_images.(layer).color = [0 0 0];
% end
end
