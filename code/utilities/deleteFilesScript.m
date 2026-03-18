% Define the main directory where the search starts
main_dir = fullfile(pwd, '..'); % Change this to your directory

% Get list of all nii files in subdirectories
files = dir(fullfile(main_dir, '**', 'PEP_WP4_EEG*_pairwise_decoding_gamma_avg_pca.mat'));

% % Create output directory if wanted
% output_dir = fullfile(pwd, '..', 'dnn_features');
% if exist(output_dir, 'dir')
%     rmdir(output_dir, 's');
% end
% 
% % Loop through found files and delete them
% for i = 1:length(files)
%     file_path = fullfile(files(i).folder, files(i).name);
%     rmdir(file_path, 's');
%     %delete(file_path);
%     fprintf('Deleted: %s\n', file_path);
% end

% Get list of all nii files in subdirectories
%files = dir(fullfile(main_dir, '**', 'typeD*'));

% Loop through found files and delete them
for i = 1:length(files)
    file_path = fullfile(files(i).folder, files(i).name);
    delete(file_path);
    fprintf('Deleted: %s\n', file_path);
end


disp('Deletion complete.');