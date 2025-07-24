function save_plot(fig, fig_name, fig_path)

% make directory if does not exists already
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

% Save the figure as a MATLAB figure (.fig) file
figFilename = fullfile(fig_path,[fig_name,'.fig']);
savefig(figFilename);

% Save the figure as a PNG image
pngFilename = fullfile(fig_path,[fig_name,'.png']);
saveas(fig, pngFilename, 'png');

end