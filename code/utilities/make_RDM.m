function [RDM_plot, mat_out, pval] = make_RDM(mat_in, cfg)

% evaluate input
if ~isfield(cfg, 'correlation_type');  cfg.correlation_type = 'spearman';end
if ~isfield(cfg, 'labels'); cfg.labels = []; end
if ~isfield(cfg, 'cell_label_style'); cfg.cell_label_style = 'none'; end
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'dissimilarity'); cfg.dissimilarity = true; end
if ~isfield(cfg, 'regressOutMean'); cfg.regressOutMean = false; end
if ~isfield(cfg, 'plot_rdm'); cfg.plot_rdm = false; end

if istable(mat_in)

    % if its a table make it an array and remove all text
    numericColumns = varfun(@isnumeric, mat_in, 'OutputFormat', 'uniform');
    mat_in =  mat_in(:, numericColumns);

    cfg.labels = mat_in.Properties.VariableNames;
    mat_in = table2array(mat_in);
end

% get group average
if cfg.regressOutMean

    % get mean
    groupMean = mean(mat_in, 2, 'omitnan');

    % loop through subjects and regress out mean
    regressedMat = zeros(size(mat_in));
    for iCol = 1:width(mat_in)
        % Design matrix: group-average and intercept
        X = [groupMean, ones(height(mat_in), 1)];
        % Perform regression
        beta = X \ mat_in(:, iCol); % Compute coefficients
        predicted = X * beta; % Predicted values based on the group average
        % Residual (indivdual column with group average regressed out)
        regressedMat(:, iCol) = mat_in(:, iCol) - predicted;
    end

    % overwrite the matrix
    mat_in = regressedMat;
end

% make correlation
[mat_out,pval] = corr(double(mat_in), 'type', cfg.correlation_type, 'Rows', 'pairwise');

% make RDM
if cfg.dissimilarity
    mat_out = 1 - mat_out;
    mat_out(eye(size(mat_out)) == 1) = 0;
    if ~isfield(cfg, 'MinColorValue'); cfg.MinColorValue = 0; end
    if ~isfield(cfg, 'MaxColorValue'); cfg.MaxColorValue = 2; end
else
    if ~isfield(cfg, 'MinColorValue'); cfg.MinColorValue = -1; end
    if ~isfield(cfg, 'MaxColorValue'); cfg.MaxColorValue = 1; end
end

% plotting
if cfg.plot_rdm
    RDM_plot = plot_RDM(mat_out, pval, cfg);
else
    RDM_plot = [];
end
end