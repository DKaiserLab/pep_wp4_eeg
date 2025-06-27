function RDM_plot = plot_RDM(mat_out, pval, cfg)

% evaluate input
if ~isfield(cfg, 'labels'); cfg.labels = [];end
if ~isfield(cfg, 'cell_label_style'); cfg.cell_label_style = 'none'; end
if ~isfield(cfg, 'colormap'); cfg.colormap = cfg.colormaps.white_zero; end
if ~isfield(cfg, 'MinColorValue'); cfg.MinColorValue = 0; end
if ~isfield(cfg, 'MaxColorValue'); cfg.MaxColorValue = 2; end
if ~isfield(cfg, 'new_figure'); cfg.new_figure = false; end

% set ticks appearance
show_all_ticks = 1;
if length(cfg.labels) > 30
    show_all_ticks = 0;
end

% generate cell labels
if strcmp(cfg.cell_label_style, 'pval')
    cell_labels = pval;
elseif strcmp(cfg.cell_label_style, '*')
    % get asterisks
    correction = 'fdr';
    cell_labels =  pval2asterisks(pval, correction);
elseif strcmp(cfg.cell_label_style, 'coef')
    cell_labels = round(mat_out, 1);
else
    cell_labels = [];
end

% make plot
if cfg.new_figure
    RDM_plot = figure;
end
% with labels
if width(mat_out) == length(cfg.labels)
    RDM_plot = deoras_heatmap(mat_out, cfg.labels, cfg.labels, cell_labels,...
        'Colormap', cfg.colormap, 'Colorbar', true, 'MinColorValue', cfg.MinColorValue, 'MaxColorValue', cfg.MaxColorValue,...
        'TickAngle', 45, 'ShowAllTicks', show_all_ticks);

    % without labels
else
    RDM_plot = deoras_heatmap(mat_out, [], [], cell_labels,...
        'Colormap', cfg.colormap, 'Colorbar', true, 'MinColorValue', cfg.MinColorValue, 'MaxColorValue', cfg.MaxColorValue,...
        'TickAngle', 45, 'ShowAllTicks', show_all_ticks);
    disp('Length of labels do not match matrix dimensions')
end

title('Correlation Plot');
drawnow

end
