function images = load_images(cfg)

% preparation
if ~isfield(cfg, 'categories');  cfg.categories = {'kitchen', 'bathroom'};end
if ~isfield(cfg, 'analysis_names');  cfg.analysis_names = {'typical', 'control', 'photos'};end
images = struct;

% get image files for drawigns
drawing_image_folder = fullfile(pwd, '..', 'drawings');
drawing_image_files = dir(fullfile(drawing_image_folder, '*gen*.png'));
drawing_image_files = struct2table(drawing_image_files);

% get image files for
photo_image_folder = fullfile(pwd, '..', 'photos');
photo_image_files = dir(fullfile(photo_image_folder, '*.png'));
photo_image_files = struct2table(photo_image_files);


for category = cfg.categories
    category = char(category);

    for analysis_name = cfg.analysis_names
        analysis_name = char(analysis_name);

        % filter for category and analysis name
        if strcmp(analysis_name, 'control')
            image_idx = (contains(drawing_image_files.name, category(1:3)) & ...
                contains(drawing_image_files.name, 'copy'));
            file_list = drawing_image_files(image_idx, :);
        elseif strcmp(analysis_name, 'typical')
            image_idx = (contains(drawing_image_files.name, category(1:3)) & ...
                ~contains(drawing_image_files.name, 'copy'));
            file_list = drawing_image_files(image_idx, :);
        elseif strcmp(analysis_name, 'photos')
            image_idx = (contains(photo_image_files.name, category(1:3)));
            file_list = photo_image_files(image_idx, :);
        end
        

        for i = 1:height(file_list)

            % get image name
            field_name = char(file_list.name(i));

            % load image
            image_path = fullfile(file_list.folder(i), field_name);
            new_image = imread(char(image_path));

            % make gray scale
            if ndims(new_image) == 3
                new_image = rgb2gray(new_image);
            end

            % Store the image data in structre
            images.(analysis_name).(category).image_data{i} = new_image;
            images.(analysis_name).(category).image_names{i} = field_name;

            % get subject numbers from images names
            images.(analysis_name).(category).sub_num_img(i) = str2num(field_name(1:3));

            % show progress
            disp(['Loaded image: ',field_name])
        end
    end
end

end

